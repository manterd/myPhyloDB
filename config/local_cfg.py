from django.conf import settings


def update():
    # Debug
    settings.DEBUG = True
    settings.ALLOWED_HOSTS = ['*']

    # Settings for django-registration
    settings.REGISTRATION_OPEN = True
    settings.ACCOUNT_ACTIVATION_DAYS = 7
    settings.REGISTRATION_EMAIL_SUBJECT_PREFIX = 'myPhyloDB Registration: '
    settings.SEND_ACTIVATION_EMAIL = False
    settings.REGISTRATION_AUTO_LOGIN = False

    # Settings for email server
    settings.EMAIL_USE_TLS = True
    settings.EMAIL_HOST = 'smtp.gmail.com'
    settings.EMAIL_PORT = 587
    settings.EMAIL_HOST_USER = 'admin@example.com'
    settings.EMAIL_HOST_PASSWORD = 'password'
