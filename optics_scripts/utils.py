from Bio import SeqIO
from Bio.Seq import Seq
from openpyxl import Workbook
from openpyxl.styles import PatternFill

STANDARD_AMINO_ACIDS = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
}

NUCLEOTIDE_IUPAC = set("ACGTURYSWKMBDHVN")
VALID_DNA_IUPAC = NUCLEOTIDE_IUPAC - {'U'}


def _letters_only(sequence):
    """Return uppercase sequence letters, preserving biological ambiguity codes."""
    return ''.join(ch for ch in str(sequence).upper() if ch.isalpha())


def _normalise_nucleotide_sequence(sequence):
    """Return an uppercase DNA sequence with RNA U converted to T."""
    return _letters_only(sequence).replace('U', 'T')


def _validate_nucleotide_sequence(sequence):
    """Raise a clear error if a sequence contains non-IUPAC nucleotide letters."""
    invalid = sorted(set(sequence) - VALID_DNA_IUPAC)
    if invalid:
        raise ValueError(f"Invalid nucleotide character(s): {', '.join(invalid)}")


def is_probable_nucleotide_sequence(sequence, min_length=651, nucleotide_fraction_threshold=0.90):
    """
    Heuristically classify an input record as nucleotide rather than protein.

    The threshold is intentionally conservative because many amino-acid one-letter
    codes overlap IUPAC nucleotide ambiguity codes. Full opsin coding sequences
    should be longer than the usual OPTICS protein-length upper bound and almost
    entirely nucleotide characters. Shorter nucleotide records can still be
    translated by explicitly setting input_seq_type='nucleotide'.
    """
    letters = _letters_only(sequence)
    if len(letters) < min_length:
        return False

    nuc_count = sum(1 for ch in letters if ch in NUCLEOTIDE_IUPAC)
    return (nuc_count / len(letters)) >= nucleotide_fraction_threshold


def reverse_complement(sequence):
    """Return the reverse complement of a DNA/RNA sequence using Biopython Seq."""
    seq = _normalise_nucleotide_sequence(sequence)
    _validate_nucleotide_sequence(seq)
    return str(Seq(seq).reverse_complement())


def translate_nucleotide_sequence(sequence, frame=1, trim_terminal_stop=True, genetic_code=1):
    """
    Translate a nucleotide sequence with Biopython's Seq.translate().

    Args:
        sequence: DNA or RNA sequence. Whitespace, digits, gaps, and FASTA line
            wrapping should already have been removed by the caller, but any
            non-letter characters are ignored defensively.
        frame: Reading frame. Use 1, 2, or 3 for the forward strand and -1, -2,
            or -3 for the reverse-complement strand.
        trim_terminal_stop: Remove one terminal stop codon from the translated
            protein if present. Internal stops are preserved as '*'.
        genetic_code: Biopython translation table identifier. Defaults to table 1
            (the standard genetic code), which is appropriate for nuclear opsin CDS input.

    Returns:
        Translated amino-acid sequence. Ambiguous codons are handled by Biopython,
        which returns 'X' for ambiguous amino-acid calls under the standard table.
    """
    try:
        frame = int(frame)
    except (TypeError, ValueError):
        raise ValueError("Translation frame must be one of 1, 2, 3, -1, -2, -3, or 'auto'.")

    if frame not in {1, 2, 3, -1, -2, -3}:
        raise ValueError("Translation frame must be one of 1, 2, 3, -1, -2, or -3.")

    seq = _normalise_nucleotide_sequence(sequence)
    _validate_nucleotide_sequence(seq)

    if frame < 0:
        seq = reverse_complement(seq)

    offset = abs(frame) - 1
    coding = seq[offset:]
    coding = coding[:len(coding) - (len(coding) % 3)]

    if not coding:
        translated = ''
    else:
        translated = str(Seq(coding).translate(table=genetic_code, to_stop=False))

    if trim_terminal_stop and translated.endswith('*'):
        translated = translated[:-1]
    return translated


def _translation_score(protein, min_protein_length=250, max_protein_length=650):
    """Score a translated protein; lower scores indicate a better frame."""
    internal_stops = protein.count('*')
    ambiguous = protein.count('X')
    standard_len = sum(1 for aa in protein if aa in STANDARD_AMINO_ACIDS)

    if min_protein_length <= standard_len <= max_protein_length:
        length_penalty = 0
    elif standard_len < min_protein_length:
        length_penalty = min_protein_length - standard_len
    else:
        length_penalty = standard_len - max_protein_length

    no_start_penalty = 0 if protein.startswith('M') else 1
    return (internal_stops, length_penalty, ambiguous, no_start_penalty, -standard_len)


def translate_nucleotide_auto_frame(sequence, min_protein_length=250, max_protein_length=650):
    """
    Translate all six reading frames and return the best candidate protein.

    Frame choice prioritizes the fewest internal stops, then a plausible opsin
    protein length, then fewer ambiguous codons, then an initiating methionine,
    and finally the longest standard-amino-acid sequence.
    """
    candidates = []
    for frame in (1, 2, 3, -1, -2, -3):
        protein = translate_nucleotide_sequence(sequence, frame=frame)
        score = _translation_score(protein, min_protein_length, max_protein_length)
        candidates.append((score, frame, protein))

    score, frame, protein = min(candidates, key=lambda item: item[0])
    return protein, frame, {
        'internal_stops': score[0],
        'length_penalty': score[1],
        'ambiguous_codons': score[2],
        'standard_aa_length': -score[4]
    }


def prepare_sequence_for_prediction(sequence, input_seq_type='auto', translation_frame='auto'):
    """
    Prepare one FASTA record for OPTICS amino-acid prediction.

    Args:
        sequence: Raw sequence body without FASTA header.
        input_seq_type: 'auto', 'protein', or 'nucleotide'. In auto mode,
            probable nucleotide records are translated and protein records are
            left unchanged.
        translation_frame: 'auto' for six-frame selection, or one of 1, 2, 3,
            -1, -2, -3.

    Returns:
        (prepared_sequence, metadata) where prepared_sequence is the amino-acid
        sequence used downstream and metadata describes any translation.
    """
    input_seq_type = (input_seq_type or 'auto').lower().strip()
    if input_seq_type not in {'auto', 'protein', 'nucleotide'}:
        raise ValueError("input_seq_type must be 'auto', 'protein', or 'nucleotide'.")

    seq = ''.join(str(sequence).upper().split())
    should_translate = input_seq_type == 'nucleotide'
    detected_type = 'protein'
    if input_seq_type == 'auto' and is_probable_nucleotide_sequence(seq):
        should_translate = True
        detected_type = 'nucleotide'
    elif input_seq_type == 'nucleotide':
        detected_type = 'nucleotide'

    metadata = {
        'detected_input_type': detected_type if input_seq_type == 'auto' else input_seq_type,
        'translated': False,
        'translation_frame': '',
        'translation_internal_stops': '',
        'translation_ambiguous_codons': '',
        'translation_method': '',
        'raw_sequence_length': len(_letters_only(seq)),
        'prepared_sequence_length': len(seq)
    }

    if not should_translate:
        return seq, metadata

    if str(translation_frame).lower().strip() == 'auto':
        protein, selected_frame, metrics = translate_nucleotide_auto_frame(seq)
    else:
        selected_frame = int(translation_frame)
        protein = translate_nucleotide_sequence(seq, frame=selected_frame)
        metrics = {
            'internal_stops': protein.count('*'),
            'ambiguous_codons': protein.count('X')
        }

    metadata.update({
        'translated': True,
        'translation_frame': selected_frame,
        'translation_internal_stops': metrics.get('internal_stops', ''),
        'translation_ambiguous_codons': metrics.get('ambiguous_codons', ''),
        'translation_method': 'Biopython Seq.translate(table=1)',
        'prepared_sequence_length': len(protein)
    })
    return protein, metadata


def extract_fasta_entries(file):
    """
    Read FASTA records with Biopython SeqIO while preserving OPTICS' historical
    return format: (names, entries), where each entry is a FASTA-formatted string.
    """
    names = []
    sequences = []

    for record in SeqIO.parse(file, 'fasta'):
        name = str(record.description).strip().replace(' ', '_').replace('\n', '')
        sequence = str(record.seq).replace('\n', '')
        names.append(name)
        sequences.append(f'>{name}\n{sequence}')

    return names, sequences


def write_to_excel(names, predictions, per_iden_list, output_filename="output.xlsx", 
                  mean_predictions=None, median_predictions=None, ci_lowers=None, 
                  ci_uppers=None, std_dev_list=None, hex_color_list=None, seq_lens_list=None):
    """
    Writes data to an Excel sheet, including bootstrap statistics and
    hexadecimal color codes, and colors the cells based on the hex codes.

    Args:
        names: List of names.
        predictions: List of predictions.
        per_iden_list: List of percentage identities.
        output_filename: Name of the output Excel file.
        mean_predictions: List of mean predictions (optional, for bootstrap).
        median_predictions: List of median predictions (optional, for bootstrap).
        ci_lowers: List of lower confidence intervals (optional, for bootstrap).
        ci_uppers: List of upper confidence intervals (optional, for bootstrap).
        std_dev_list: List of standard deviations (optional, for bootstrap).
        hex_color_list: List of hexadecimal color codes.
        seq_lens_list: List of sequence lengths
        """

    wb = Workbook()
    ws = wb.active

    
    if mean_predictions == None:
        ws.append(['Names', 'Single_Prediction', '%Identity_Nearest_VPOD_Sequence', 'Sequence_Length','Lmax_Hex_Color'])
        for i in range(len(names)):
            # Because openpyxel is picky about hex-codes we need to remove the '#' symbol for it to accept it as a fill color.
            hex_color = hex_color_list[i].replace('#','') 
            ws.append([names[i], predictions[i], per_iden_list[i], seq_lens_list[i], hex_color_list[i]])
            ws.cell(row=i+2, column=5).fill = PatternFill(start_color=hex_color, 
                                                        end_color=hex_color, 
                                                        fill_type="solid")
    else:
        ws.append(['Names', 'Single_Prediction', 'Prediction_Means', 'Prediction_Medians',
                    'Prediction_Lower_Bounds', 'Prediction_Upper_Bounds', 'Std_Deviation', 
                    '%Identity_Nearest_VPOD_Sequence', 'Sequence_Lengths','Lmax_Hex_Color'])
        for i in range(len(names)):
            # Because openpyxel is picky about hex-codes we need to remove the '#' symbol for it to accept it as a fill color.
            hex_color = hex_color_list[i].replace('#','') 
            ws.append([names[i], predictions[i], mean_predictions[i], median_predictions[i],
                        ci_lowers[i], ci_uppers[i], std_dev_list[i], per_iden_list[i], seq_lens_list[i], hex_color_list[i]])
            ws.cell(row=i+2, column=10).fill = PatternFill(start_color=hex_color, 
                                                        end_color=hex_color, 
                                                        fill_type="solid")
    wb.save(output_filename)
