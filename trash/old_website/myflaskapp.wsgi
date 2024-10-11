import sys
import os

sys.path.insert(0, '/var/www/lmax_pred')  # Adjust the path accordingly

activate_this = '~/home/.local/share/virtualenvs/optics/bin/scripts/activate_this.py'
exec(open(activate_this).read(), dict(__file__=activate_env)) 

from app import app as application  # Assuming your Flask app variable is 'app'
