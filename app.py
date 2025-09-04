import sys

from flask import (Flask, jsonify, render_template, request,
                   send_from_directory, session)
from functions import WebApp


app = Flask(__name__)
app.secret_key = "super_secret_key"

# Initialize: Application
webapp = WebApp()

# Figure storage
figures = {}


from flask_wtf.csrf import CSRFProtect, CSRFError

csrf = CSRFProtect(app)


@app.route('/run', methods=['POST'])
def run():
    data = request.get_json() # Get the JSON body
    message = data.get('message') # Extract the message

    # Call method
    data = webapp.pressButton(message)

    return jsonify(data)

    # enzymeName = request.form.get('enzymeName')
    # fileExp = request.files.get('fileExp')



@app.route('/evalDNA', methods=['POST'])
def evalDNA():
    print('Form')
    data = {}
    data['enzymeName'] = request.form.get('enzymeName')
    data['seq5Prime'] = request.form.get('seq5Prime')
    data['seq3Prime'] = request.form.get('seq3Prime')
    data['seqLength'] = request.form.get('seqLength')
    data['minPhred'] = request.form.get('minPhred')

    data['fileExp'] = request.files.get('fileExp')
    data['fileExpRev'] = request.files.get('fileExpRev')  # add this
    data['fileBg'] = request.files.get('fileBg')
    data['fileBgRev'] = request.files.get('fileBgRev')  # add this

    session['figures'] = webapp.evalDNA(data)

    # Return JSON directly instead of redirecting
    return jsonify(figures)


@app.route('/data/figures/<filename>')
def getFigure(filename):
    return send_from_directory('data/figures', filename)


@app.route('/checkFigures')
@csrf.exempt
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


if __name__ == '__main__':
    app.run(debug=True)
