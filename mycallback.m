function [res] = mycallback(obj,event_obj)
% obj          handle to the object that has been clicked on
% event_obj    handle to event data object (empty in this release)
% res [output] logical flag to determine whether the rotate
%              operation should take place or the 'ButtonDownFcn'
%              property of the object should take precedence
%turn pan off
assignin('base','obj',obj)
event_obj.Axes
			rotate3d(event_obj.Axes,'off')

disp('app')

						ggg= event_obj.Axes.View
			%sprintf('Current view az %d el %d',ggg)
			azel=ggg;
			%assignin('base','gg',event_obj.Axes.View)
			%azel = reshape(ggg,[2 1]);
			azel(1) = mod(azel(1),360)-180;
			ttmp=azel2phitheta(reshape(azel,[2,1]),false);

			%sprintf('Current view phi %d theta %d',ttmp)
			%ttmp=azel2phitheta([azel(1); azel(2)],false);
ttmp
			disp('get phi theta')
			app.PhiEditField.Value = mod(ttmp(1)-270,360);
			app.ThetaEditField.Value = ttmp(2);

			rotate3d(event_obj.Axes,'off')

			%Calculate(app)
			res=true;
end