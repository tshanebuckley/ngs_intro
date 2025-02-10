from os import listdir
from os.path import isfile, join

def get_messages():
    # string returned form the request
    # default is for if no messages were found
    unavailable = 'No messages are available.'

    # attempt to find the available messages
    messages_path = './messages'
    text_files = [f[:-4] for f in listdir(messages_path) if (isfile(join(messages_path, f)) and f.endswith('.txt'))]
    ret_str = "\n".join(text_files)

    # if there were no message files
    if not ret_str:
        return unavailable

    # return the request
    return ret_str

def get_message(message):
    # string returned form the request
    # default is for if no message was found
    ret_str = f'No message found for "{message}".'
    
    # attempt to return the message requested    
    path = f'./messages/{message}.txt'
    if isfile(path):
        with open(path, "r") as file:
            lines = file.readlines()
            if len(lines) > 0:
                ret_str = lines[0]

    # return the request
    return ret_str