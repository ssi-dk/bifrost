. activate bifrost
uwsgi -s /tmp/run_checker.sock --manage-script-name --mount /=run_checker:application --http :8051