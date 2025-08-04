from flask import Flask, jsonify, render_template, request
from functions import WebApp
import json


app = Flask(__name__)

# Initialize: Application
webapp = WebApp()


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
    data = request.get_json() # Get the input form
    webapp.evalDNA(data)

    return jsonify({'status': 'ok'})


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
