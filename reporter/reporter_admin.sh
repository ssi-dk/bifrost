. activate bifrost
uwsgi -s /tmp/reporter_admin.sock --manage-script-name --mount /=reporter_admin:server --http :8052 --ini reporter_admin.ini