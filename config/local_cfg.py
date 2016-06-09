from django.conf import settings


# User-defined settings for multi-threading (integer values only)
usr_num_threads = 3     # number of analyses that can be run simultaneously
usr_numcore = 4    # number of processors per analysis


# Django settings
def update():
    settings.DEBUG = True
    settings.ALLOWED_HOSTS = ['*']
    settings.TEMPLATE_DEBUG = False
    settings.SESSION_EXPIRE_AT_BROWSER_CLOSE = True
    settings.USE_L10N = True
    settings.USE_I18N = True
    settings.LANGUAGE_CODE = 'en-us'
    settings.USE_TZ = True
    settings.TIME_ZONE = 'America/Denver'
    settings.REGISTRATION_OPEN = True

