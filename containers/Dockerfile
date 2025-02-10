# Use a base image
FROM docker.io/python:3.9

# Set the working directory
WORKDIR /app

# Make the directory where we expect our messages to be mounted
RUN mkdir /app/messages

# Copy the application code
COPY ./hello /app

# Install dependencies (flask and click)
RUN pip install Flask
RUN pip install click

# Expose a port (if necessary)
EXPOSE 80

# Define the command to run when the container starts
CMD ["python", "app.py"]