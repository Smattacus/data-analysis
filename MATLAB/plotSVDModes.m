function plotSVDModes(filename)
%Quick function to plot the SVD modes of a given file.
%
% plotSVDModes(filename)
%
%
[X V S U] = genSVDFromArray(filename);
%Generate the polar components for the topos plot
p = (0:2 * pi / 8 : (2 * pi - 2 * pi / 8))';
p = [p; p(1)];
%Solid line for reference at r = 1.
s = (0:2 * pi / 1000: 2 * pi);
t = ones(size(s,1), size(s,2));
%Plot the different modes - sinusoidal pattern is pretty apparent.
if sum(find([U(:,1); U(1,1)] > 1)) == 0 
    figure(1); hold; 
    polar(s, t, '-black');
    polar(p, 1 + [U(:,1); U(1,1)],'-square'); title('First SVD Mode');
else
    figure(1); 
    polar(p, 1 + [U(:,1); U(1,1)],'-square'); title('First SVD Mode');
    hold; polar(s, t, '-black');hold;
end
figure(2); polar(p, 1 + [U(:,2); U(1,2)],'-square'); title('Second SVD Mode');
hold; polar(s, t, '-black');hold;
figure(3); polar(p, 1 + [U(:,3); U(1,3)],'-square'); title('Third SVD Mode');
hold; polar(s, t, '-black');hold;
figure(4); polar(p, 1 + [U(:,4); U(1,4)],'-square'); title('Fourth SVD Mode');
hold; polar(s, t, '-black');hold;
figure(5); polar(p, 1 + [U(:,5); U(1,5)],'-square'); title('Fifth SVD Mode');
hold; polar(s, t, '-black');hold;
figure(6); polar(p, 1 + [U(:,6); U(1,6)],'-square'); title('Sixth SVD Mode');
hold; polar(s, t, '-black');hold;
figure(7); polar(p, 1 + [U(:,7); U(1,7)],'-square'); title('Seventh SVD Mode');
hold; polar(s, t, '-black');hold;
figure(8); polar(p, 1 + [U(:,8); U(1,8)],'-square'); title('Eighth SVD Mode');
hold; polar(s, t, '-black');hold;
