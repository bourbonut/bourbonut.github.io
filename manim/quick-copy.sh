PRJROOT=".."
MEDIA=$PRJROOT/manim/media
POST=$PRJROOT/posts/involute-spur-gear

# force
force=$MEDIA/images/force/Gear2Gear_ManimCE_v0.17.2.png
cp $force $POST/force.png

# involute
static_involute=$MEDIA/images/involute/Static_ManimCE_v0.17.2.png
video_involute=$MEDIA/videos/involute/720p30/InvoluteFunction.mp4
cp $static_involute $POST/static_involute.png
cp $video_involute $POST

# movement
gearonrack=$MEDIA/videos/movement/720p30/GearOnRack.mp4
rackongear=$MEDIA/videos/movement/720p30/RackOnGear.mp4
closerackongear=$MEDIA/videos/movement/720p30/CloseRackOnGear.mp4
cp $gearonrack $POST
cp $rackongear $POST
cp $closerackongear $POST

# rack_gear
gearwithinterference=$MEDIA/images/rack_gear/GearWithInterference_ManimCE_v0.17.2.png
gearwithoutinterference=$MEDIA/images/rack_gear/GearWithoutInterference_ManimCE_v0.17.2.png
rackandgear=$MEDIA/images/rack_gear/RackAndGear_ManimCE_v0.17.2.png
rack=$MEDIA/images/rack_gear/Rack_ManimCE_v0.17.2.png
cp $gearwithinterference $POST/gear_with_interference.png
cp $gearwithoutinterference $POST/gear_without_interference.png
cp $rackandgear $POST/rack_and_gear.png
cp $rack $POST/rack.png
