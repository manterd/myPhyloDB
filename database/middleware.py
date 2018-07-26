from django.conf import settings
from django.core.cache import cache
from importlib import import_module


class UserRestrictMiddleware(object):
    def process_request(self, request):
        """
        Checks if different session exists for user and deletes it.
        """
        if request.user.is_authenticated():
            thiscache = cache
            cache_timeout = 86400
            cache_key = "user_pk_%s_restrict" % request.user.pk
            cache_value = thiscache.get(cache_key)

            if cache_value is not None:
                if request.session.session_key != cache_value:
                    engine = import_module(settings.SESSION_ENGINE)
                    session = engine.SessionStore(session_key=cache_value)
                    session.delete()
                    thiscache.set(cache_key, request.session.session_key,
                              cache_timeout)
            else:
                thiscache.set(cache_key, request.session.session_key, cache_timeout)