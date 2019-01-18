. activate bifrost
gunicorn --bind '0.0.0.0:8051' run_checker:server