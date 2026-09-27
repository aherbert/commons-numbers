% Licensed to the Apache Software Foundation (ASF) under one or more
% contributor license agreements.  See the NOTICE file distributed with
% this work for additional information regarding copyright ownership.
% The ASF licenses this file to You under the Apache License, Version 2.0
% (the "License"); you may not use this file except in compliance with
% the License.  You may obtain a copy of the License at
%
%     https://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,\
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Matlab hurwitzZeta function requires the Symbolic Math Toolbox

% Read a file of (s, a) values, evaluate zeta(s, a), and write to a result file.
% Change variable precision arithmetic (VPA) precision using function digits(n).
function hzeta(filename)
  M = readmatrix(filename,'CommentStyle','#');
  if size(M, 2) < 2
    error('Input file must contain at least two columns: s and a.');
  end
  [dir, name] = fileparts(filename);
  fileID = fopen(fullfile(dir, [name, '.csv']), 'w');
  fprintf(fileID, "# Evaluated using Matlab %s\n", version);
  fprintf(fileID, "# hurwitzZeta(s, a) using VPA digits=%d\n", digits);
  for i=1:size(M,1)
    s = M(i,1);
    a = M(i,2);
    z = vpa(hurwitzZeta(sym(s, 'f'), sym(a, 'f')));
    fprintf(fileID, "%.17g, %.17g, %s\n", s, a, z);
  end
  fclose(fileID);
end
