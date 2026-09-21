# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,\
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import argparse
import sys
import os
from mpmath import mp, zeta, __version__ as version
parser = argparse.ArgumentParser(description="Program to compute zeta(s).")
parser.add_argument("data", type=str, help="s values")
parser.add_argument("--negate", default=False, action=argparse.BooleanOptionalAction, help="negate s")
args = parser.parse_args()
mp.dps = 36
mp.pretty = True

# Read the header and duplicate it
with open(sys.argv[0]) as f:
  for line in f:
    if line[0] != '#':
      break
    print(line, end='')

print()
print(f'# Evaluated using {os.path.basename(sys.argv[0])}')
print(f'# mpmath ({version}) zeta')
print()

with open(args.data) as f:
  for line in f:
    if line[0] == '#':
      print(line, end='')
      continue
    s = float(line)
    if args.negate:
      s = -s
    print(f'{s}, {zeta(s)}')
