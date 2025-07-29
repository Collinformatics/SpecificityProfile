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
    const form = document.getElementById("dnaForm");
    const formData = new FormData(form);

    fetch('/evalDNA', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        console.log("Response from server:", data);
        // You can now update the page dynamically with JS
    })
    .catch(error => {
        console.error("Error:", error);
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

