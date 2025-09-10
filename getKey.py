import secrets

def genKey(app):
    app.config['SECRET_KEY'] = secrets.token_hex(nbytes=32)  # required for CSRF
    print(f'Key: {app.config['SECRET_KEY']}\n'
          f'Len: {len(app.config['SECRET_KEY'])}\n')