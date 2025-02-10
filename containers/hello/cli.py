import click
from util import get_messages, get_message

@click.group()
def group():
    pass

@click.command()
def hello():
    click.echo('Hello, World!')

@click.command()
def messages():
    click.echo(get_messages())

@click.command()
@click.argument('msg')
def message(msg):
    click.echo(get_message(msg))

group.add_command(hello)
group.add_command(messages)
group.add_command(message)

if __name__ == '__main__':
    group()