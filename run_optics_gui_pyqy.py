import sys
import os
import math
import random
import traceback
import warnings

# PyQt6 Imports
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, 
                             QComboBox, QCheckBox, QGroupBox, QScrollArea, 
                             QTextEdit, QFileDialog, QMessageBox, QDialog, 
                             QSizePolicy, QFrame)
from PyQt6.QtCore import Qt, QTimer, QThread, pyqtSignal, QObject, QSize, QRectF
from PyQt6.QtGui import QPainter, QColor, QPen, QBrush, QFont, QIcon, QPalette

warnings.filterwarnings("ignore")

# --- Backend Imports & Mocking ---
# We keep these safe so the GUI runs even if the backend files are missing
try:
    from optics_predictions import run_optics_predictions
except ImportError:
    run_optics_predictions = None

try:
    from optics_shap import generate_shap_explanation
except ImportError:
    generate_shap_explanation = None

# --- Custom Logging Stream ---
class StreamRedirector(QObject):
    """Redirects stdout/stderr to a PyQt Signal so it can be shown in the GUI safely."""
    text_written = pyqtSignal(str)

    def write(self, text):
        self.text_written.emit(str(text))

    def flush(self):
        pass

# --- Eye Animation Loading Screen ---
class HumanEyeLoadingScreen(QDialog):
    """
    A modal loading screen with an animated eye drawn using QPainter.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.WindowType.Dialog | Qt.WindowType.FramelessWindowHint)
        self.setModal(True)
        self.setFixedSize(400, 300)
        self.setStyleSheet("background-color: #1c1c1c; border: 1px solid #333;")

        # Animation State
        self.gaze_current_x = 0.0
        self.gaze_target_x = 0.0
        self.gaze_current_y = 0.0
        self.gaze_target_y = 0.0
        
        self.gaze_positions = [
            (-40, -25), (45, 20), (0, 0), (30, 25), (-20, -20), (0, 0), 
            (-45, 0), (45, 0), (0, 30), (0, -30), (0, 0)
        ]
        self.gaze_position_index = 0
        self.gaze_timer_counter = 100
        self.ellipsis_count = 0

        # Timers
        self.anim_timer = QTimer(self)
        self.anim_timer.timeout.connect(self.update_animation)
        
        self.text_timer = QTimer(self)
        self.text_timer.timeout.connect(self.update_text)

    def showEvent(self, event):
        super().showEvent(event)
        self.start_animation()

    def start_animation(self):
        self.anim_timer.start(30) # ~30fps
        self.text_timer.start(400)

    def stop_animation(self):
        self.anim_timer.stop()
        self.text_timer.stop()
        self.accept()

    def update_text(self):
        self.ellipsis_count = (self.ellipsis_count + 1) % 4
        self.update() # Trigger repaint to draw text

    def update_animation(self):
        self.gaze_timer_counter -= 1
        
        # Pick new target
        if self.gaze_timer_counter <= 0:
            self.gaze_position_index = (self.gaze_position_index + 1) % len(self.gaze_positions)
            self.gaze_target_x, self.gaze_target_y = self.gaze_positions[self.gaze_position_index]
            
            if self.gaze_target_x == 0 and self.gaze_target_y == 0:
                self.gaze_timer_counter = random.randint(100, 150)
            else:
                self.gaze_timer_counter = random.randint(40, 80)

        # Smooth easing
        self.gaze_current_x += (self.gaze_target_x - self.gaze_current_x) * 0.1
        self.gaze_current_y += (self.gaze_target_y - self.gaze_current_y) * 0.1
        
        self.update() # Trigger paintEvent

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        # Center coordinates
        center_x = self.width() / 2
        center_y = self.height() / 2 - 20 # Shift up slightly

        # 1. Draw Eye Sclera (White part)
        painter.setBrush(QBrush(QColor("#EAEAEA")))
        painter.setPen(Qt.PenStyle.NoPen)
        # Outer dark rim
        painter.setBrush(QBrush(QColor("#111111")))
        painter.drawEllipse(QRectF(center_x - 110, center_y - 80, 220, 160))
        # White inner
        painter.setBrush(QBrush(QColor("#EAEAEA")))
        painter.drawEllipse(QRectF(center_x - 100, center_y - 70, 200, 140))

        # 2. Draw Iris
        iris_x = center_x + self.gaze_current_x
        iris_y = center_y + self.gaze_current_y
        iris_radius = 35
        painter.setBrush(QBrush(QColor("#5DADE2")))
        painter.drawEllipse(QRectF(iris_x - iris_radius, iris_y - iris_radius, iris_radius*2, iris_radius*2))

        # 3. Draw Pupil
        pupil_radius = 15
        painter.setBrush(QBrush(QColor("black")))
        painter.drawEllipse(QRectF(iris_x - pupil_radius, iris_y - pupil_radius, pupil_radius*2, pupil_radius*2))

        # 4. Highlight
        highlight_x = iris_x + 10
        highlight_y = iris_y - 10
        painter.setBrush(QBrush(QColor("white")))
        painter.drawEllipse(QRectF(highlight_x - 5, highlight_y - 5, 8, 8))

        # 5. Draw Text
        dots = "." * self.ellipsis_count
        text = f"Processing{dots}"
        font = QFont("Century Gothic", 14, QFont.Weight.Bold)
        painter.setFont(font)
        painter.setPen(QColor("white"))
        
        text_rect = QRectF(0, center_y + 90, self.width(), 40)
        painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, text)


# --- Worker Thread for Long Running Tasks ---
class WorkerThread(QThread):
    finished_signal = pyqtSignal(bool, str) # success, message

    def __init__(self, mode, params):
        super().__init__()
        self.mode = mode
        self.params = params

    def run(self):
        try:
            if self.mode == 'predictions':
                if run_optics_predictions is None:
                    raise ImportError("run_optics_predictions module not found.")
                
                # Unpack params
                pred_df, output_path = run_optics_predictions(
                    input_sequence=self.params['input_sequence'],
                    pred_dir=self.params['pred_dir'],
                    output=self.params['output'],
                    model=self.params['model'],
                    encoding_method=self.params['encoding_method'],
                    blastp=self.params['blastp'],
                    iden_report=self.params['iden_report'],
                    refseq=self.params['refseq'],
                    reffile=self.params['reffile'],
                    bootstrap=self.params['bootstrap'],
                    visualize_bootstrap=self.params['visualize_bootstrap'],
                    bootstrap_viz_file=self.params['bootstrap_viz_file'],
                    save_as=self.params['save_as'],
                    full_spectrum_xaxis=self.params['full_spectrum_xaxis'],
                    model_version=self.params['model_version'],
                    tolerate_non_standard_aa=self.params['tolerate_non_standard_aa'],
                    tolerate_incomplete_seqs=self.params['tolerate_incomplete_seqs']
                )
                msg = f"Results are in: {os.path.dirname(output_path)}" if output_path else "Run complete."
                self.finished_signal.emit(True, f"OPTICS predictions completed successfully!\n{msg}")

            elif self.mode == 'shap':
                if generate_shap_explanation is None:
                    raise ImportError("generate_shap_explanation module not found.")

                generate_shap_explanation(
                    input_file=self.params['input_file'], 
                    pred_dir=self.params['pred_dir'], 
                    output=self.params['output'], 
                    save_as=self.params['save_as'],
                    model=self.params['model'], 
                    encoding_method=self.params['encoding_method'], 
                    model_version=self.params['model_version'], 
                    cmd_line="GUI_Execution", 
                    mode=self.params['shap_mode'], 
                    n_positions=self.params['n_positions'], 
                    use_reference_sites=self.params['use_reference_sites']
                )
                self.finished_signal.emit(True, f"SHAP analysis completed successfully!\nResults are in: {self.params['pred_dir']}")

        except Exception as e:
            tb = traceback.format_exc()
            print(f"\n--- ERROR ---\n{tb}") # This goes to the GUI log via stdout redirect
            self.finished_signal.emit(False, str(e))


# --- Mode Selector Dialog ---
class ModeSelectorDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Select OPTICS Mode")
        self.setFixedSize(500, 400)
        self.selected_mode = None
        
        layout = QVBoxLayout(self)
        layout.setSpacing(20)
        layout.setContentsMargins(40, 40, 40, 40)

        # Logo/Title
        title_lbl = QLabel("Welcome to OPTICS")
        title_lbl.setStyleSheet("font-family: 'Century Gothic'; font-size: 24px; font-weight: bold;")
        title_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title_lbl)

        sub_lbl = QLabel("Select your analysis pipeline:")
        sub_lbl.setStyleSheet("font-family: 'Century Gothic'; font-size: 14px;")
        sub_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(sub_lbl)

        # Buttons
        btn_style = """
            QPushButton {
                font-family: 'Century Gothic'; 
                font-size: 16px; 
                padding: 15px;
                background-color: #3E3E3E;
                color: white;
                border-radius: 5px;
            }
            QPushButton:hover { background-color: #505050; }
        """

        pred_btn = QPushButton("Standard Predictions\n(Spectral Tuning)")
        pred_btn.setStyleSheet(btn_style)
        pred_btn.clicked.connect(lambda: self.select_mode('predictions'))
        layout.addWidget(pred_btn)

        shap_btn = QPushButton("SHAP Interpretation\n(Feature Importance)")
        shap_btn.setStyleSheet(btn_style)
        shap_btn.clicked.connect(lambda: self.select_mode('shap'))
        layout.addWidget(shap_btn)

    def select_mode(self, mode):
        self.selected_mode = mode
        self.accept()


# --- Main Window ---
class OpticsMainWindow(QMainWindow):
    def __init__(self, mode):
        super().__init__()
        self.mode = mode
        self.init_ui()
        self.init_logging()

    def init_ui(self):
        title_suffix = "Predictions" if self.mode == 'predictions' else "SHAP Analysis"
        self.setWindowTitle(f"OPTICS: {title_suffix}")
        self.resize(1000, 950)
        
        # Central Widget & Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # --- Top Bar ---
        top_bar = QHBoxLayout()
        header_lbl = QLabel(f"OPTICS: {title_suffix}")
        header_lbl.setStyleSheet("font-family: 'Century Gothic'; font-size: 20px; font-weight: bold;")
        top_bar.addWidget(header_lbl)
        
        top_bar.addStretch()
        
        self.theme_btn = QPushButton("Toggle Theme") # Placeholder, logic handled by Palette mostly
        self.theme_btn.clicked.connect(self.toggle_theme)
        top_bar.addWidget(self.theme_btn)
        main_layout.addLayout(top_bar)

        # --- Scroll Area ---
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        content_widget = QWidget()
        self.form_layout = QVBoxLayout(content_widget)
        self.form_layout.setSpacing(15)
        scroll.setWidget(content_widget)
        main_layout.addWidget(scroll)

        # --- Common Inputs ---
        common_group = QGroupBox("General Options")
        common_layout = QVBoxLayout() # Use VBox for rows of HBox
        common_group.setLayout(common_layout)
        
        # Helper to create file inputs
        def create_file_input(label_text, btn_cmd):
            row = QHBoxLayout()
            lbl = QLabel(label_text)
            lbl.setFixedWidth(180)
            entry = QLineEdit()
            btn = QPushButton("Browse...")
            btn.clicked.connect(btn_cmd)
            row.addWidget(lbl)
            row.addWidget(entry)
            row.addWidget(btn)
            common_layout.addLayout(row)
            return entry

        # Helper for combos
        def create_combo_input(label_text, choices):
            row = QHBoxLayout()
            lbl = QLabel(label_text)
            lbl.setFixedWidth(180)
            combo = QComboBox()
            combo.addItems(choices)
            row.addWidget(lbl)
            row.addWidget(combo)
            common_layout.addLayout(row)
            return combo

        self.input_file_edit = create_file_input("Input Sequence/FASTA:", self.browse_input)
        self.output_dir_edit = create_file_input("Output Directory:", self.browse_output)
        
        # Output Prefix
        row_prefix = QHBoxLayout()
        lbl_prefix = QLabel("Output Filename Prefix:")
        lbl_prefix.setFixedWidth(180)
        self.prefix_edit = QLineEdit()
        row_prefix.addWidget(lbl_prefix)
        row_prefix.addWidget(self.prefix_edit)
        common_layout.addLayout(row_prefix)

        # Combos
        self.version_choices = ['vpod_1.3']
        self.model_choices = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 
                              'wildtype-vert', 'type-one', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 
                              'wildtype-vert-mnm', 'wildtype-mut']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.viz_ftyp_choices = ['svg', 'png', 'pdf']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']

        self.version_combo = create_combo_input("Model Version:", self.version_choices)
        self.model_combo = create_combo_input("Model:", self.model_choices)
        self.encoding_combo = create_combo_input("Encoding Method:", self.encoding_choices)
        self.encoding_combo.setCurrentIndex(1) # aa_prop default

        self.form_layout.addWidget(common_group)

        # --- Mode Specific Widgets ---
        if self.mode == 'predictions':
            self.setup_predictions_ui()
        else:
            self.setup_shap_ui()

        # --- Run Button ---
        btn_text = "Run OPTICS Predictions" if self.mode == 'predictions' else "Run SHAP Analysis"
        self.run_btn = QPushButton(btn_text)
        self.run_btn.setMinimumHeight(50)
        self.run_btn.setStyleSheet("font-size: 16px; font-weight: bold; margin-top: 10px;")
        self.run_btn.clicked.connect(self.start_process)
        self.form_layout.addWidget(self.run_btn)

        # --- Log Output ---
        log_group = QGroupBox("Output Log")
        log_layout = QVBoxLayout()
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setMinimumHeight(200)
        self.log_output.setStyleSheet("font-family: Consolas, monospace; font-size: 10pt;")
        log_layout.addWidget(self.log_output)
        log_group.setLayout(log_layout)
        self.form_layout.addWidget(log_group)

    def setup_predictions_ui(self):
        # Toggles
        self.tol_non_std_cb = QCheckBox("Tolerate non-standard Amino Acids")
        self.tol_non_std_cb.setChecked(True)
        self.form_layout.addWidget(self.tol_non_std_cb)

        self.tol_incomp_cb = QCheckBox("Tolerate Incomplete Sequences (<250 or >650 AA)")
        self.form_layout.addWidget(self.tol_incomp_cb)

        # BLASTp Options
        self.blast_group = QGroupBox("BLASTp Options")
        self.blast_group.setCheckable(True)
        self.blast_group.setChecked(True)
        blast_layout = QFormLayout()
        
        self.blast_report_edit = QLineEdit("blastp_report.txt")
        blast_layout.addRow("BLASTp Report Filename:", self.blast_report_edit)

        self.refseq_combo = QComboBox()
        self.refseq_combo.addItems(self.refseq_choices)
        self.refseq_combo.currentTextChanged.connect(self.toggle_custom_ref)
        blast_layout.addRow("Reference Sequence:", self.refseq_combo)

        ref_file_layout = QHBoxLayout()
        self.custom_ref_edit = QLineEdit()
        self.custom_ref_btn = QPushButton("Browse...")
        self.custom_ref_btn.clicked.connect(self.browse_custom_ref)
        ref_file_layout.addWidget(self.custom_ref_edit)
        ref_file_layout.addWidget(self.custom_ref_btn)
        blast_layout.addRow("Custom Reference File:", ref_file_layout)
        
        self.blast_group.setLayout(blast_layout)
        self.form_layout.addWidget(self.blast_group)
        
        # Trigger initial state
        self.toggle_custom_ref(self.refseq_combo.currentText())

        # Bootstrap Options
        self.boot_group = QGroupBox("Bootstrap Options")
        self.boot_group.setCheckable(True)
        self.boot_group.setChecked(True)
        boot_layout = QVBoxLayout()

        self.viz_boot_cb = QCheckBox("Visualize Bootstrap Predictions")
        self.viz_boot_cb.setChecked(True)
        self.viz_boot_cb.toggled.connect(self.toggle_boot_viz)
        boot_layout.addWidget(self.viz_boot_cb)

        # Viz sub-options
        self.viz_container = QWidget()
        viz_layout = QFormLayout(self.viz_container)
        self.boot_viz_edit = QLineEdit("bootstrap_viz")
        viz_layout.addRow("Bootstrap Viz Filename:", self.boot_viz_edit)
        
        self.viz_type_combo = QComboBox()
        self.viz_type_combo.addItems(self.viz_ftyp_choices)
        viz_layout.addRow("Viz Filetype:", self.viz_type_combo)

        self.full_xaxis_cb = QCheckBox("Enable Full-Spectrum X-axis (300-650nm)")
        viz_layout.addRow(self.full_xaxis_cb)
        
        boot_layout.addWidget(self.viz_container)
        self.boot_group.setLayout(boot_layout)
        self.form_layout.addWidget(self.boot_group)

    def setup_shap_ui(self):
        shap_group = QGroupBox("SHAP Interpretation Options")
        shap_layout = QFormLayout()

        self.shap_mode_combo = QComboBox()
        self.shap_mode_combo.addItems(['both', 'comparison', 'single'])
        shap_layout.addRow("Analysis Mode:", self.shap_mode_combo)

        self.n_pos_edit = QLineEdit("10")
        shap_layout.addRow("Top N Features:", self.n_pos_edit)

        self.use_ref_cb = QCheckBox("Use Reference Numbering")
        shap_layout.addRow(self.use_ref_cb)

        self.shap_save_combo = QComboBox()
        self.shap_save_combo.addItems(self.viz_ftyp_choices)
        shap_layout.addRow("Save Visualizations As:", self.shap_save_combo)

        shap_group.setLayout(shap_layout)
        self.form_layout.addWidget(shap_group)

    def init_logging(self):
        # Redirect stdout/stderr
        self.redirector = StreamRedirector()
        self.redirector.text_written.connect(self.append_log)
        sys.stdout = self.redirector
        sys.stderr = self.redirector

    def append_log(self, text):
        self.log_output.moveCursor(self.log_output.textCursor().MoveOperation.End)
        self.log_output.insertPlainText(text)
        self.log_output.ensureCursorVisible()

    # --- Interaction Slots ---

    def toggle_theme(self):
        # This is a basic switch. Ideally, use pyqtdarktheme or a full stylesheet.
        # We implemented a default Dark Fusion palette below in main().
        # Reverting to light is tricky without storing the original palette, 
        # but for this snippet we'll just toggle a flag (functionality omitted for brevity 
        # as standard Qt palettes are global).
        QMessageBox.information(self, "Theme", "Theme toggling requires restarting or advanced palette management in PyQt.")

    def browse_input(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Select Input FASTA", "", "FASTA (*.fasta *.fa *.fna);;Text (*.txt);;All (*.*)")
        if fname: self.input_file_edit.setText(fname)

    def browse_output(self):
        dname = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dname: self.output_dir_edit.setText(dname)

    def browse_custom_ref(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Select Custom Reference", "", "FASTA (*.fasta *.fa);;All (*.*)")
        if fname: self.custom_ref_edit.setText(fname)

    def toggle_custom_ref(self, text):
        enabled = (text == "custom" and self.blast_group.isChecked())
        self.custom_ref_edit.setEnabled(enabled)
        self.custom_ref_btn.setEnabled(enabled)

    def toggle_boot_viz(self, checked):
        self.viz_container.setEnabled(checked)

    # --- Processing ---

    def start_process(self):
        # Validation
        if not self.input_file_edit.text():
            QMessageBox.critical(self, "Input Error", "Please specify an input file.")
            return

        # Prepare Params
        params = {
            'input_file': self.input_file_edit.text(), # For SHAP
            'input_sequence': self.input_file_edit.text(), # For Preds
            'pred_dir': self.output_dir_edit.text() or os.path.join(os.getcwd(), 'prediction_outputs'),
            'output': self.prefix_edit.text() or "optics_results",
            'model_version': self.version_combo.currentText(),
            'model': self.model_combo.currentText(),
            'encoding_method': self.encoding_combo.currentText(),
        }

        if self.mode == 'predictions':
            params.update({
                'tolerate_non_standard_aa': self.tol_non_std_cb.isChecked(),
                'tolerate_incomplete_seqs': self.tol_incomp_cb.isChecked(),
                'blastp': self.blast_group.isChecked(),
                'iden_report': self.blast_report_edit.text(),
                'refseq': self.refseq_combo.currentText(),
                'reffile': self.custom_ref_edit.text() if self.refseq_combo.currentText() == "custom" else None,
                'bootstrap': self.boot_group.isChecked(),
                'visualize_bootstrap': self.viz_boot_cb.isChecked(),
                'bootstrap_viz_file': self.boot_viz_edit.text(),
                'save_as': self.viz_type_combo.currentText(),
                'full_spectrum_xaxis': self.full_xaxis_cb.isChecked()
            })
        else:
            try:
                n_pos = int(self.n_pos_edit.text())
            except ValueError:
                n_pos = 10
            
            params.update({
                'shap_mode': self.shap_mode_combo.currentText(),
                'n_positions': n_pos,
                'use_reference_sites': self.use_ref_cb.isChecked(),
                'save_as': self.shap_save_combo.currentText()
            })

        # UI State
        self.run_btn.setEnabled(False)
        self.log_output.clear()
        print(f"Starting {self.mode.upper()} analysis...")

        # Loading Screen
        self.loading_screen = HumanEyeLoadingScreen(self)
        self.loading_screen.show()

        # Thread
        self.worker = WorkerThread(self.mode, params)
        self.worker.finished_signal.connect(self.on_process_finished)
        self.worker.start()

    def on_process_finished(self, success, message):
        self.loading_screen.stop_animation()
        self.loading_screen.close()
        self.run_btn.setEnabled(True)
        
        if success:
            QMessageBox.information(self, "Success", message)
        else:
            QMessageBox.critical(self, "Error", f"An error occurred:\n{message}")


def set_dark_theme(app):
    """Manually apply a Fusion Dark Theme."""
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.ColorRole.WindowText, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ColorRole.ToolTipBase, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.ToolTipText, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ColorRole.ButtonText, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.BrightText, Qt.GlobalColor.red)
    palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.ColorRole.HighlightedText, Qt.GlobalColor.black)
    app.setPalette(palette)


def main():
    app = QApplication(sys.argv)
    
    # Set global font
    font = QFont("Century Gothic", 10)
    app.setFont(font)
    
    set_dark_theme(app)

    # 1. Selector
    selector = ModeSelectorDialog()
    if selector.exec() == QDialog.DialogCode.Accepted:
        # 2. Main Window
        mode = selector.selected_mode
        window = OpticsMainWindow(mode)
        window.show()
        sys.exit(app.exec())
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()