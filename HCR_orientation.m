function HCR_angle=HCR_orientation(w)
count=0;
for angle=-90:90
    count=count+1;
    data=squeeze(w(:,:,t))';
    data=imrotate(data,angle,'bilinear','crop');
    roll_flux=size(data,2)*mean(data,2);
    for i=1:size(data,1)
        line=squeeze(data(i,:));
        if roll_flux(i)>0
            MFR(i)=roll_flux(i)/sum(line(line>0));
        elseif roll_flux(i)<0
            MFR(i)=roll_flux(i)/sum(line(line<0));
        end
    end
    mean_MFR(count)=mean(MFR);
end
[maxval,maxloc]=max(mean_MFR);
HCR_angle=-90+maxloc;