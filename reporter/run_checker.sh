uwsgi -s /tmp/reporter.sock --manage-script-name --mount /=run_checker:application --http :8051
