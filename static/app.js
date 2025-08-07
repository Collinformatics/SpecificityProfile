// Fix AA
function updateFixedAA() {
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const container = document.getElementById('fixedAAContainer');
    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];

    container.innerHTML = '';
    container.style.display = 'flex';
    container.style.flexWrap = 'wrap'; // Wrap to next line if needed
    container.style.gap = '12px'; // spacing between each checkbox
    container.paddingLeft = '5px';

    for (let i = 1; i <= seqLength; i++) {
        const wrapper = document.createElement('div');
        wrapper.style.display = 'flex';
        wrapper.style.flexDirection = 'column';
        wrapper.style.flex = '0 0 60px';

        const label = document.createElement('label');
        label.style.display = 'flex';
        label.style.alignItems = 'center';
        label.style.color = '#FFF'; // '#FA8128'
        label.style.flex = '0 0 30px';

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = `filterPos`;
        checkbox.value = `R${i}`;

        const text = document.createTextNode(`R${i}`);
        label.appendChild(checkbox);
        label.appendChild(text);
        wrapper.appendChild(label); // Append label to wrapper

        const aaGroup = document.createElement('div');
        aaGroup.style.display = 'none';
        aaGroup.style.flexWrap = 'wrap';
        aaGroup.style.gap = '0px';
        aaGroup.style.marginLeft = '10px';
        aaGroup.style.marginTop = '0px';

        aminoAcids.forEach(aa => {
            const aaLabel = document.createElement('label');
            aaLabel.style.display = 'flex';
            aaLabel.style.alignItems = 'center';
            aaLabel.style.color = 'white';
            aaLabel.style.fontSize = '12px';

            const aaCheckbox = document.createElement('input');
            aaCheckbox.type = 'checkbox';
            aaCheckbox.name = `fixR${i}`;
            aaCheckbox.value = aa;

            const aaText = document.createTextNode(aa);
            aaLabel.appendChild(aaCheckbox);
            aaLabel.appendChild(aaText);
            aaGroup.appendChild(aaLabel);
        });

        checkbox.addEventListener('change', function () {
            aaGroup.style.display = checkbox.checked ? 'flex' : 'none';
        });

        wrapper.appendChild(aaGroup);      // Append AA group to wrapper
        container.appendChild(wrapper);    // Append full wrapper to container
    }
}

// Create listener to inspect html input after the form has been loaded and parsed
document.addEventListener('DOMContentLoaded', function() {
    updateFixedAA();

    // Inspect File Upload: Experimental data
    document.getElementById('fileExp').addEventListener('change', function(event) {
    const file = event.target.files[0];
    if (!file) return;

    const allowed = ['.fasta', '.fa', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fastq.gz', '.fq.gz'];
    const filename = file.name.toLowerCase();

    const valid = allowed.some(ext => filename.endsWith(ext));
    if (!valid) {
        alert(`ERROR:\nYou can only upload files with these extensions:
        ${allowed.join(' ')}`);
        event.target.value = ''; // Clear file input
    }
    });

    // Inspect File Upload: Background data
    document.getElementById('fileBg').addEventListener('change', function(event) {
    const file = event.target.files[0];
    if (!file) return;

    const allowed = ['.fasta', '.fa', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fastq.gz', '.fq.gz'];
    const filename = file.name.toLowerCase();

    const valid = allowed.some(ext => filename.endsWith(ext));
    if (!valid) {
        alert(`ERROR:\nYou can only upload files with these extensions:
        ${allowed.join(' ')}`);
        event.target.value = ''; // Clear file input
    }
    });

//    // Inspect File Upload:
//    document.getElementById('').addEventListener('input', function() {
//        const value = this.value;
//        console.log('Input: ', value)
//
//        if (!isInteger && value !== '') {
//            alert('Please enter a valid input.');
//            this.value = ''; // Optionally clear the input
//        }
//    });
//});


//// Form: Process DNA
//function buttonProcessDNA_UpdateWhenDoneBuilding() {
//    const form = document.getElementById("dnaForm");
//    const formData = new FormData(form); // This keeps files intact!
//
//    fetch('/evalDNA', {
//        method: 'POST',
//        body: formData // send as multipart/form-data
//    })
//    .then(response => response.json())
//    .then(data => {
//        console.log("Response from server:", data);
//    })
//    .catch(error => {
//        console.error("Error:", error);
//    });
//
//    // app.py
//    @app.route('/evalDNA', methods=['POST'])
//    def evalDNA():
//        enzymeName = request.form.get('enzymeName')
//        fileExp = request.files.get('fileExp')
//        fileBg = request.files.get('fileBg')
//        # other form fields the same way
//
});


// Define button function
function buttonProcessDNA() {
    const form = document.getElementById("dnaForm");
    const formData = new FormData(form);
    const json = {}; // Dont send files as a JSON
    const selectedFixPositions = [];

    // Process the input form
    for (const [key, value] of formData.entries()) {
    if (key === 'filterPos') {
        selectedFixPositions.push(value);  // e.g., ['R2']
        }


        if (json[key]) {
        // When you have more that one value or a key, put the values in a list
            if (!Array.isArray(json[key])) {
                json[key] = [json[key]]; // Convert to array
            }
            json[key].push(value);  // Push another value into the list
        } else {
            json[key] = value;
        }
        }

        // Clean out fixR* keys not selected
        Object.keys(json).forEach(key => {
            if (key.startsWith('fix') && !selectedFixPositions.includes(key.replace('fix', ''))) {
                delete json[key];
        }
    });


    console.log('Input Form:', json);
        fetch('/evalDNA', {
        method: 'POST',
        body: formData  // Send the actual FormData object, not a JSON
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


//
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
