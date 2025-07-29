// Define button function
function clickButton() {
    const message = 'your data';  // Pass data to app.py

    // POST request to app.py
    fetch('/run', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json' // Send JSON (Optional)
        },
        body: JSON.stringify({ message: message }) // Send data to app.py
    })
    .then(response => response.text())
    .then(data => {
        console.log('Server Response:\n', data);
    })
    .catch(error => {
        console.error('ERROR: ', error);
    });
}

function buttonProcessDNA() {
    // Get the input value
    const enzymeName = document.getElementById("enzymeName").value;

    // POST request to app.py
    fetch('/evalDNA', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json' // Send JSON (Optional)
        },
        body: JSON.stringify({ enzymeName: enzymeName }) // Send data to app.py
    })
    .then(response => response.text())
    .then(data => {
        console.log('Server Response:\n', data);
    })
    .catch(error => {
        console.error('ERROR: ', error);
    });
}


function pageHome() {
    window.location.href = "/"
}

function pageProcessDNA() {
    window.location.href = "/processDNA";
}

function pageFilterAA() {
    window.location.href = "/filterAminoAcids";
}

function pageFilterMotif() {
    window.location.href = "/filterMotif";
}

