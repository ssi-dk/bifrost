. activate bifrost
gunicorn --bind '0.0.0.0:8052' reporter_admin:server