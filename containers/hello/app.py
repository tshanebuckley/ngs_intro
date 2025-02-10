from flask import Flask
from util import get_messages, get_message

app = Flask(__name__)

@app.route('/')
def hello():
    return 'Hello, World!'

@app.route('/messages')
def my_messages():
    return get_messages()

@app.route('/messages/<message>')
def my_message(message):
    return get_message(message)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=80)