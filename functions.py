class WebApp:
    def __init__(self):
        self.buttonState = False
        self.messages = []

    def pressButton(self, message):
        print(f'Received data: {message}')

        return {'key': 'Returned data'}


    def evalDNA(self, data):
        print(f'Received Data:\n{data}\n\n')

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}
        enzymeName = data['enzymeName']
        fileExp = data['fileExp']  # This is a FileStorage object


        # Read the file contents (as bytes)
        contents = fileExp.read().decode('utf-8')  # Decode to string if text-based

        print(f'enzymeName: {enzymeName}')
        print(f'File contents:\n{contents[:200]}')

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}

    def filterAA(self):
        print('Processing Substrates')

        return {'AA': 'VEHTVALKQNR'}

    def filterMotif(self):
        print('Filtering Substrates')

        return {'Motif': 'TVALK'}
