from flask_cache import Cache
from flask_debugtoolbar import DebugToolbarExtension
from flask_assets import Environment
from flask_wtf.csrf import CSRFProtect

# Setup flask cache
cache = Cache()

# init flask assets
assets_env = Environment()

# Toolbar for debugging
debug_toolbar = DebugToolbarExtension()

# To enable CSRF protection globally for a Flask app
csrf = CSRFProtect()
 
