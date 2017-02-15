from django.conf import settings


# number of analyses that can be run simultaneously
usr_num_threads = 3


# Django settings
def update():
    settings.DEBUG = False
    settings.ALLOWED_HOSTS = ['*']
    settings.TEMPLATE_DEBUG = False
    settings.SESSION_EXPIRE_AT_BROWSER_CLOSE = True
    settings.USE_L10N = True
    settings.USE_I18N = True
    settings.LANGUAGE_CODE = 'en-us'
    settings.USE_TZ = True
    settings.TIME_ZONE = 'America/Denver'
