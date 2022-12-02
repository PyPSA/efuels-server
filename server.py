## Copyright 2018-2020 Tom Brown

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU Affero General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## License and more information at:
## https://github.com/PyPSA/whobs-server




from flask import Flask, render_template

app = Flask(__name__)


#defaults to only listen to GET and HEAD
@app.route('/')
def root():
    return render_template('index.html')


if __name__ == '__main__':
    app.run(port='5002')
