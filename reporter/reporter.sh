. activate bifrost
gunicorn --bind '0.0.0.0:8050' reporter:server