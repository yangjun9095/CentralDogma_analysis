
%% Plot for mRNA boundary pos, width
%Plot for xpos for different Dm
plot(Dm/60,Result(:,2))

%Plot for width for different Dm
%plot(Dm,Result(:,2)

% title('Position of the Boundary for different Dm','fontsize',25)
% xlabel('mRNA Diffusion Constant(um^2/sec)','fontsize',20)
% ylabel('Boundary Position in AP axis','fontsize',20)

title('Width of the Boundary for different Dm','fontsize',25)
xlabel('mRNA Diffusion Constant(um^2/sec)','fontsize',20)
ylabel('Boundary Width in AP axis','fontsize',20)

%% Plot for Protein boundary pos, width
%Plot for xpos for different Dp
plot(Dp/60,ResultP(:,1))

%Plot for width for different Dp
%plot(Dp/60,ResultP(:,2))

title('Position of the Boundary for different Dp','fontsize',25)
xlabel('Protein Diffusion Constant(um^2/sec)','fontsize',20)
ylabel('Boundary Position in AP axis','fontsize',20)
ylim([0 1])
% title('Width of the Boundary for different Dp','fontsize',25)
% xlabel('Protein Diffusion Constant(um^2/sec)','fontsize',20)
% ylabel('Boundary Width in AP axis','fontsize',20)