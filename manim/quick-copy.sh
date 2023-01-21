# manim -pqm involute.py -a
# manim -pqm movement.py -a
# manim -pqm force.py -a
# manim -pqm rack_gear.py -a

PRJROOT=".."
MEDIA=$PRJROOT/manim/media
SPUR=$PRJROOT/posts/involute-spur-gear
SNN=$PRJROOT/posts/introduction-spiking-nn

# force
force=$MEDIA/images/force/Gear2Gear_ManimCE_v0.17.2.png
cp $force $SPUR/force.png

# involute
static_involute=$MEDIA/images/involute/Static_ManimCE_v0.17.2.png
video_involute=$MEDIA/videos/involute/720p30/InvoluteFunction.mp4
cp $static_involute $SPUR/static_involute.png
cp $video_involute $SPUR

# movement
gearonrack=$MEDIA/videos/movement/720p30/GearOnRack.mp4
rackongear=$MEDIA/videos/movement/720p30/RackOnGear.mp4
closerackongear=$MEDIA/videos/movement/720p30/CloseRackOnGear.mp4
cp $gearonrack $SPUR
cp $rackongear $SPUR
cp $closerackongear $SPUR

# rack_gear
gearwithinterference=$MEDIA/images/rack_gear/GearWithInterference_ManimCE_v0.17.2.png
gearwithoutinterference=$MEDIA/images/rack_gear/GearWithoutInterference_ManimCE_v0.17.2.png
rackandgear=$MEDIA/images/rack_gear/RackAndGear_ManimCE_v0.17.2.png
rack=$MEDIA/images/rack_gear/Rack_ManimCE_v0.17.2.png
cp $gearwithinterference $SPUR/gear_with_interference.png
cp $gearwithoutinterference $SPUR/gear_without_interference.png
cp $rackandgear $SPUR/rack_and_gear.png
cp $rack $SPUR/rack.png

# signal
impulse=$MEDIA/videos/signal/720p30/Impulse.mp4
cp $impulse $SNN
