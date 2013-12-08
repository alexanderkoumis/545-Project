% a quick hack to generate the SL link parameters from Nao specifications

% read the original Nao inertial parameters in the coordinates given by Aldebaran
fp = fopen('NaoOriginalInertialParameters.txt','r');

parms = {};
c = 0;

while 1,
  [name,count] = fscanf(fp,'%s',1);
  [s,count] = fscanf(fp,'%f  %f %f %f  %f %f %f %f %f %f %f %f %f',13);
  if count < 13,
    break;
  else
    c=c+1;
    parms(c).name      = name;
    parms(c).m         = s(1)/1000;
    parms(c).cm        = s(2:4)/1000;
    parms(c).inertia   = reshape(s(5:13),3,3)/1000/1000/1000;
  end
end

fclose(fp);

% some fudging for numerial stability
parms(7).inertia(1,1) = 0.001;
parms(6).inertia(1,1) = 0.001;
parms(5).inertia(1,1) = 0.005;

% give the AAA joint too much inertia to stabilize the entire anke
parms(9).inertia(1,1) = parms(9).inertia(1,1)*100;
parms(9).inertia(2,2) = parms(9).inertia(2,2)*100;
parms(9).inertia(3,3) = parms(9).inertia(3,3)*100;
parms(10).inertia(3,3) = parms(10).inertia(3,3)*100;
parms(14).inertia(1,1) = parms(14).inertia(1,1)*10;
parms(14).inertia(2,2) = parms(14).inertia(2,2)*10;
parms(14).inertia(3,3) = parms(14).inertia(3,3)*10;
parms(13).inertia(1,1) = parms(13).inertia(1,1)*100;
parms(12).inertia(2,2) = parms(12).inertia(2,2)*10;

% the joint ordering from SL
c = 0;

c = c+1; linkparms(c).name = 'BASE';

c = c+1; linkparms(c).name = 'R_SFE';
c = c+1; linkparms(c).name = 'R_SAA';
c = c+1; linkparms(c).name = 'R_HR';
c = c+1; linkparms(c).name = 'R_EB';
c = c+1; linkparms(c).name = 'R_WR';
c = c+1; linkparms(c).name = 'R_FING';

c = c+1; linkparms(c).name = 'L_SFE';
c = c+1; linkparms(c).name = 'L_SAA';
c = c+1; linkparms(c).name = 'L_HR';
c = c+1; linkparms(c).name = 'L_EB';
c = c+1; linkparms(c).name = 'L_WR';
c = c+1; linkparms(c).name = 'L_FING';

c = c+1; linkparms(c).name = 'R_FB';
c = c+1; linkparms(c).name = 'R_HFE';
c = c+1; linkparms(c).name = 'R_HAA';
c = c+1; linkparms(c).name = 'R_KFE';
c = c+1; linkparms(c).name = 'R_AFE';
c = c+1; linkparms(c).name = 'R_AAA';

c = c+1; linkparms(c).name = 'L_FB';
c = c+1; linkparms(c).name = 'L_HFE';
c = c+1; linkparms(c).name = 'L_HAA';
c = c+1; linkparms(c).name = 'L_KFE';
c = c+1; linkparms(c).name = 'L_AFE';
c = c+1; linkparms(c).name = 'L_AAA';

c = c+1; linkparms(c).name = 'B_HR';
c = c+1; linkparms(c).name = 'B_HN';

% set all rotation matrices for the joints to eye by default
for i=1:length(linkparms),
  linkparms(i).R = eye(3);
end

% set R for special joints
linkparms(1).R  = rotmat_axis([0;0;1],-pi/2);  % the base coordinates

linkparms(14).R = rotmat_axis([1;0;0],-pi/4);  % the forebend coordinates
linkparms(20).R = rotmat_axis([1;0;0],-pi/4);

M = [1 -1 1; -1 1 -1; 1 -1 1]; % matrix to mirror an inertia in the Nao y-axis
m = [1; -1; 1];

% search through the orignal parms and copy data over
for i=1:length(linkparms)
  for j=1:length(parms)
    if strcmp(parms(j).name,linkparms(i).name), % these are right side parameters
      disp(linkparms(i).name);
      linkparms(i).inertia = linkparms(i).R'*parms(j).inertia*linkparms(i).R;
      linkparms(i).mcm = parms(j).m*linkparms(i).R'*parms(j).cm;
      linkparms(i).m   = parms(j).m;
    elseif strcmp(parms(j).name(2:end),linkparms(i).name(2:end)), % infer left side parameters
      disp(linkparms(i).name);
      linkparms(i).inertia = linkparms(i).R'*(parms(j).inertia.*M)*linkparms(i).R;
      linkparms(i).mcm = parms(j).m*linkparms(i).R'*(parms(j).cm.*m);
      linkparms(i).m   = parms(j).m;
    end
  end
end

% the finger enerial we have to make up
linkparms(7).m = 0.2;
linkparms(7).mcm = [0.03;0;0];
linkparms(7).inertia = linkparms(7).m*...
  (eye(3)*(linkparms(7).mcm'*linkparms(7).mcm)-linkparms(7).mcm*linkparms(7).mcm')+...
  eye(3)*0.001;
linkparms(7).mcm = linkparms(7).m*linkparms(7).mcm;

linkparms(13) = linkparms(7);
linkparms(13).name = 'L_FING';

% print all link parameters in SL format
for i=1:length(linkparms),
  disp(sprintf('%s   %f    %f %f %f   %f %f %f %f %f %f  %f %f %f %f',...
    linkparms(i).name,...
    linkparms(i).m,...
    linkparms(i).mcm(1),...
    linkparms(i).mcm(2),...
    linkparms(i).mcm(3),...
    linkparms(i).inertia(1,1),...
    linkparms(i).inertia(1,2),...
    linkparms(i).inertia(1,3),...
    linkparms(i).inertia(2,2),...
    linkparms(i).inertia(2,3),...
    linkparms(i).inertia(3,3),...
    0.05,0,0,0));
end



