from django.conf import settings


# number of analyses that can be run simultaneously
usr_num_threads = 3


# Django settings
def update():
    settings.DEBUG = True
    settings.ALLOWED_HOSTS = ['*']
    settings.TEMPLATE_DEBUG = True
    settings.SESSION_EXPIRE_AT_BROWSER_CLOSE = True
    settings.USE_L10N = True
    settings.USE_I18N = True
    settings.LANGUAGE_CODE = 'en-us'
    settings.USE_TZ = True
    settings.TIME_ZONE = 'America/Denver'

    ### Settings for django-registration-redux
    settings.REGISTRATION_OPEN = True

    '''
    ### Do not change this code (it is currently disabled)
    ### Settings for two-step authentification
    settings.ACCOUNT_ACTIVATION_DAYS = 7
    settings.REGISTRATION_EMAIL_SUBJECT_PREFIX = 'myPhyloDB Registration: '
    settings.SEND_ACTIVATION_EMAIL = True
    settings.REGISTRATION_AUTO_LOGIN = False

    ### Settings for email server
    settings.EMAIL_USE_TLS = True
    settings.EMAIL_HOST = 'smtp.gmail.com'
    settings.EMAIL_PORT = 587
    settings.EMAIL_HOST_USER = ''
    settings.EMAIL_HOST_PASSWORD = ''
    '''
