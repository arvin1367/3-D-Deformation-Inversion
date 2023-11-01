function mkpt(data,date)
load(data);
imagesc(ve_cm_v);
colormap jet
[x,y] = getpts;
x1=round(y(1));y1=round(x(1));
x=x1;y=y1;
me=zeros(1,size(ve_cm,3));
stde=zeros(1,size(ve_cm,3));
mn=zeros(1,size(ve_cm,3));
stdn=zeros(1,size(ve_cm,3));
mu=zeros(1,size(ve_cm,3));
stdu=zeros(1,size(ve_cm,3));
for i=3:size(ve_cm,3)
    te=[ve_cm(x-2,y-2,i) ve_cm(x-2,y-1,i) ve_cm(x-2,y,i) ve_cm(x-2,y+1,i) ve_cm(x-2,y+2,i)
        ve_cm(x-1,y-2,i) ve_cm(x-1,y-1,i) ve_cm(x-1,y,i) ve_cm(x-1,y+1,i) ve_cm(x-1,y+2,i)
        ve_cm(x,y-2,i) ve_cm(x,y-1,i) ve_cm(x,y,i) ve_cm(x,y+1,i) ve_cm(x,y+2,i)
        ve_cm(x+1,y-2,i) ve_cm(x+1,y-1,i) ve_cm(x+1,y,i) ve_cm(x+1,y+1,i) ve_cm(x+1,y+2,i)
        ve_cm(x+2,y-2,i) ve_cm(x+2,y-1,i) ve_cm(x+2,y,i) ve_cm(x+2,y+1,i) ve_cm(x+2,y+2,i)];
        
    tn=[vn_cm(x-2,y-2,i) vn_cm(x-2,y-1,i) vn_cm(x-2,y,i) vn_cm(x-2,y+1,i) vn_cm(x-2,y+2,i)
        vn_cm(x-1,y-2,i) vn_cm(x-1,y-1,i) vn_cm(x-1,y,i) vn_cm(x-1,y+1,i) vn_cm(x-1,y+2,i)
        vn_cm(x,y-2,i) vn_cm(x,y-1,i) vn_cm(x,y,i) vn_cm(x,y+1,i) vn_cm(x,y+2,i)
        vn_cm(x+1,y-2,i) vn_cm(x+1,y-1,i) vn_cm(x+1,y,i) vn_cm(x+1,y+1,i) vn_cm(x+1,y+2,i)
        vn_cm(x+2,y-2,i) vn_cm(x+2,y-1,i) vn_cm(x+2,y,i) vn_cm(x+2,y+1,i) vn_cm(x+2,y+2,i)];
    
    tu=[vu_cm(x-2,y-2,i) vu_cm(x-2,y-1,i) vu_cm(x-2,y,i) vu_cm(x-2,y+1,i) vu_cm(x-2,y+2,i)
        vu_cm(x-1,y-2,i) vu_cm(x-1,y-1,i) vu_cm(x-1,y,i) vu_cm(x-1,y+1,i) vu_cm(x-1,y+2,i)
        vu_cm(x,y-2,i) vu_cm(x,y-1,i) vu_cm(x,y,i) vu_cm(x,y+1,i) vu_cm(x,y+2,i)
        vu_cm(x+1,y-2,i) vu_cm(x+1,y-1,i) vu_cm(x+1,y,i) vu_cm(x+1,y+1,i) vu_cm(x+1,y+2,i)
        vu_cm(x+2,y-2,i) vu_cm(x+2,y-1,i) vu_cm(x+2,y,i) vu_cm(x+2,y+1,i) vu_cm(x+2,y+2,i)];

    me(1,i)=mean(te,'all');mn(1,i)=mean(tn,'all');mu(1,i)=mean(tu,'all');
    stde(1,i)=std(te,1,'all');stdn(1,i)=std(tn,1,'all');stdu(1,i)=std(tu,1,'all');
           
    clear te tu tn
end

errorbar(me,2*stde,'-o','LineWidth',1','MarkerSize',3,'CapSize',4,'Color',[197/255 86/255 89/255]);
hold on
errorbar(mn,2*stdn,'-o','LineWidth',1','MarkerSize',3,'CapSize',4,'Color',[84/255 172/255 117/255]);
hold on
errorbar(mu,2*stdu,'-o','LineWidth',1','MarkerSize',3,'CapSize',4,'Color',[117/255 114/255 181/255]);
hold on
set(gcf,'unit','centimeters','position',[10 10 20 10])
set(gca,'linewidth',1.5,'Fontname','Airl','fontsize',12);
ylabel('Displacement (m)');
xlabel('Date');
legend({'East','North','Vertical'},'FontSize',12,'Location','northwest')
legend('boxoff')
xlim([-5,i+5]);
set(gca,'XTick',[-1,round(i/3),round(2*i/3),i], ...
    'xticklabel',{datestr(date(1),'yyyy/dd/mm'),datestr(date(round(i/3)),'yyyy/dd/mm'),datestr(date(round(2*i/3)),'yyyy/dd/mm'),datestr(date(i),'yyyy/dd/mm')})
nstr=[num2str(x),'_',num2str(x),'_pt'];
print(nstr,'-dpng','-r600');

end
