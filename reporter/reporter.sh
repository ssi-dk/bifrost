. activate bifrost
uwsgi -s /tmp/reporter.sock --processes 3 --manage-script-name --mount /=reporter:server --http :8050