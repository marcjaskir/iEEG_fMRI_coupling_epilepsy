def create_pwd_file(username, password, fname=None):
  if fname is None:
    fname = "ieeglogin.bin"
  with open(fname, 'wb') as f:
    f.write(password.encode())
  files.download(fname)
  print("-- -- IEEG password file saved -- --\n")