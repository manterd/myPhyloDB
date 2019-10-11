# this file is separate from other utils files for the sake of circular imports:
#   in theory this is used by all files, including ones utils_df imports

flag = True    # only need to switch this one variable for the whole code to go into verbose mode
# TODO 1.4 Admin account auto-debug, switch back to request-linked debug flags.
# TODO 1.4 Error log saves original request for admin to rerun later (with admin auto debug)

def debug(*msg):    # *msg allows any number of arguments to be passed in, since we just want to print everything
    # print msg if debug flag is set to true
    if flag:
        print msg
    # would be somewhat convenient to check memory usage here
