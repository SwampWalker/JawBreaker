import sys

if len(sys.argv) != 3:
  print("usage: " + sys.argv[0] + " eos.yp0 eos.t00")
  print("\tMerges T=0 equation of state with yp=0 equation of state.")
  sys.exit(0)

# Open the files.
yp0File = open(sys.argv[1])
yp0lines = yp0File.readlines()
# Only want the T=0 block of the yp=0 eos.
yp0t0lines = yp0lines[4:115]

t0File = open(sys.argv[2])
t0lines = t0File.readlines()

# Print the header.
sys.stdout.write(t0lines.pop(0))
sys.stdout.write(t0lines.pop(0))
sys.stdout.write(t0lines.pop(0))
sys.stdout.write(t0lines.pop(0))

for line in yp0t0lines:
  sys.stdout.write(line)

for line in t0lines:
  sys.stdout.write(line)
