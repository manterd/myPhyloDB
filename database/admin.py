from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User

from database.models import UserProfile, Project, Reference


class UserProfileInline(admin.StackedInline):
    model = UserProfile


class UserProfileAdmin(UserAdmin):
    inlines = (UserProfileInline, )

admin.site.unregister(User)
admin.site.register(User, UserProfileAdmin)
admin.site.register(Project)
admin.site.register(Reference)
