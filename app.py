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
        pg1='This website will allow you to process enzyme profiling data.',
        pg2='Upload a FASTQ or FASTA file to translate the DNA sequences. The prevalence '
            'of amino acids in your datasets will be evaluated. Both the counts matrix, '
            'and protein sequences can then be downloaded and used for subsequent '
            'analysis',
        pg3='By comparing the probability of finding a given residue in the '
            'experimental (Exp) and background (Bg) datasets the Enrichment Score (ES) '
            'of each residue can be evaluated.',
        eqt1='ES =  log<sub>2</sub>(prob<sub>AA Exp</sub> / prob<sub>AA Bg</sub>)',
        pg4='Words',
        eqt2='∆S = S<sub>Max</sub> - S<sub>Shannon</sub> = log<sub>2</sub>(20) - '
             '∑(-prob<sub>AA</sub> * log<sub>2</sub>(prob<sub>AA</sub>))'
    )

@app.route('/processDNA')
def processDNA():
    return render_template(
        'processDNA.html',
        pg1='This program translates DNA sequences. The produced datasets can be '
            'download and used for subsequent analysis.'
    )

@app.route('/filterAA')
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
