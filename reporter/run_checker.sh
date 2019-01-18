. activate bifrost
uwsgi -s /tmp/run_checker.sock --manage-script-name --mount /=run_checker:server --http :8051