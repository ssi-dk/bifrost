. activate bifrost
export REPORTER_ADMIN=True
uwsgi -s /tmp/reporter.sock --processes 3 --manage-script-name --mount /=reporter:server --http :8050