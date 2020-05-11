import reporter
import dash_auth
import keys

reporter.ADMIN = True

app = reporter.app

dash_auth.BasicAuth(
    app,
    keys.USERNAME_PASSWORD
)

server = app.server