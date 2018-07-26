from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User as Users

from database.models import UserProfile, Project, Reference, Sample, Soil, UserDefined


class UserProfileInline(admin.StackedInline):
    model = UserProfile


class UserProfileAdmin(UserAdmin):
    inlines = (UserProfileInline, )


admin.site.unregister(Users)
admin.site.register(Users, UserProfileAdmin)
admin.site.register(Project)
admin.site.register(Reference)
admin.site.register(Sample)
admin.site.register(Soil)
admin.site.register(UserDefined)

