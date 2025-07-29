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
    return render_template('home.html',)

@app.route('/processDNA')
def processDNA():
    return render_template(
        'processDNA.html',
        pg1='This program translates DNA sequences. The produced datasets can be '
            'download and used for subsequent analysis.'
    )

@app.route('/filterAminoAcids')
def filterAA():
    return render_template(
        'filterAA.html',
        pg1='This program allows you to visualise your profiling data.<br><br>'
            'You can filter the data by locking one or more amino acids at a given '
            'position in the sequence.'
    )

@app.route('/filterMotif')
def filterMotif():
    return render_template('filterMotif.html',
        pg1='Apply an automated, entropy based filter to your data.'
    )


if __name__ == '__main__':
    app.run(debug=True)
