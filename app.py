from flask import Flask, jsonify, render_template, request
from functions import WebApp


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


@app.route('/')
def home():
    return render_template(
        'home.html',
        pg1='This website will allow you to process FASTQ files and extract protein '
            'sequences from the data. Amino acids in the extracted sequences are counted '
            'and by comparing the "experimental" and "background" datasets the '
            'enrichment of each residue can be evaluated.'
    )

@app.route('/processDNA')
def processDNA():
    return render_template('processDNA.html')

@app.route('/filterAA')
def filterAA():
    return render_template('filterAA.html')

@app.route('/filterMotif')
def filterMotif():
    return render_template('filterMotif.html')


if __name__ == '__main__':
    app.run(debug=True)
