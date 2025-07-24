class WebApp:
    def __init__(self):
        self.buttonState = False
        self.messages = []

    def pressButton(self, message):
        print(f'Received data: {message}')

        return {'key': 'Returned data'}

    def evalDNA(self):
        print('Processing DNA')

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}

    def filterAA(self):
        print('Processing Substrates')

        return {'AA': 'VEHTVALKQNR'}

    def filterMotif(self):
        print('Filtering Substrates')

        return {'Motif': 'TVALK'}
