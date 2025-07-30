function updateFixedAA() {
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const container = document.getElementById('fixedAAContainer');

    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];


    container.innerHTML = '';
    container.style.display = 'flex';
    container.style.flexWrap = 'wrap'; // Wrap to next line if needed
    container.style.gap = '20px'; // spacing between each checkbox
    for (let i = 1; i <= seqLength; i++) {
        const label = document.createElement('label');
        label.style.display = 'flex';
        label.style.alignItems = 'center';
        label.style.gap = '5px';
        label.style.color = '#FA8128'; // orange label

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = `fixR${i}`;
        checkbox.value = `R${i}`;

        const text = document.createTextNode(`R${i}`);

        label.appendChild(checkbox);
        label.appendChild(text);
        container.appendChild(label);
    }
}
document.addEventListener('DOMContentLoaded', function() {
    updateFixedAA();
});


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

