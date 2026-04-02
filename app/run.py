#!/usr/bin/env python

import os

from iso import create_app
from werkzeug.contrib.fixers import ProxyFix

app = create_app('iso.settings.ProdConfig')
app.wsgi_app = ProxyFix(app.wsgi_app)


if __name__ == "__main__":
    app.run(host='0.0.0.0')
