from flask import (Flask, jsonify, render_template, request,
                   send_from_directory, session)
from functions import WebApp
import secrets
import sys



app = Flask(__name__)
app.config['SECRET_KEY'] = secrets.token_hex(nbytes=32) # required for CSRF
print(f'Key: {app.config['SECRET_KEY']}\nLen: {len(app.config['SECRET_KEY'])}\n')

# Initialize: Application
webapp = WebApp()

# Figure storage
figures = {}


from flask_wtf.csrf import CSRFProtect, CSRFError, generate_csrf
csrf = CSRFProtect(app)
@app.route("/getToken")
def getToken():
    token = generate_csrf()
    print(f"CSRF Token: {token}")  # will print in server logs
    return f"Token: {token}"



def parseForm():
    # Parse the form
    data = {}
    for key, value in request.form.items():
        data[key] = value

    for key, value in request.files.items():
        if value:
            data[key] = value

    return data




@app.route('/run', methods=['POST'])
def run():
    data = request.get_json() # Get the JSON body
    message = data.get('message') # Extract the message

    # Call method
    data = webapp.pressButton(message)

    return jsonify(data)

    # enzymeName = request.form.get('enzymeName')
    # fileExp = request.files.get('fileExp')



@app.route('/evalFormDNA', methods=['POST'])
def evalDNA():
    # Parse the form
    data = parseForm()

    # Process the data
    webapp.evalDNA(data)

    return render_template('results.html')


@app.route('/data/figures/<filename>')
def getFigure(filename):
    return send_from_directory('data/figures', filename)


@app.route('/checkFigures')
def checkFigures():
    figures = session.get('figures', {})
    return jsonify(figures)


@app.route('/')
def home():
    # return render_template('home.html')
    return render_template('processDNA.html')


@app.route('/processDNA')
def processDNA():
    return render_template('processDNA.html')

@app.route('/filterAminoAcids')
def filterAA():
    return render_template('filterAA.html')

@app.route('/filterMotif')
def filterMotif():
    return render_template('filterMotif.html')


# Run the app
if __name__ == '__main__':
    app.run(debug=True)
