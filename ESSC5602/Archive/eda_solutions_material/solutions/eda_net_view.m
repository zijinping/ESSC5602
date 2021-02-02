
% view function makes a picture of a neural net
% you must set the figure and clear it before the
% call, e.g.
% figure(2); clf;
function null = eda_net_view(N,w,b)

Nmax = max(N);
Lmax = length(N);

axis([ 0 (Lmax+1) 0 (Nmax+1) ]);
hold on;

scale = zeros(Lmax,1);
for x = 1:Lmax
    scale(x,1) =  Nmax/N(x);
end


for x = 1:Lmax
    for y = 1:N(x)
        %plot(x,y,'o','MarkerSize',50);
        rectangle('Position',[x-.2,y*scale(x,1)-.4,.4,.8],'LineWidth',2);
        text(x,y*scale(x,1),strcat('b =  ', num2str( b(y,x))), 'HorizontalAlignment' , 'center');
        if x ~= Lmax
            for y2 = 1:N(x+1)
                if w( y2 , y , x+1) ~= 0
                    plot([x+.2 x+.8],[y*scale(x,1) y2*scale(x+1,1)],'Color',[0.7 0.7 0.7]);
                    text(x+.5,(y*scale(x,1)+y2*scale(x+1,1))/2,strcat('w =  ', num2str( w( y2 , y , x+1))), 'HorizontalAlignment' , 'center');
                end    
            end
        end
    end
end
end