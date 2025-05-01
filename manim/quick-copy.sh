# manim -pqm involute.py -a
# manim -pqm movement.py -a
# manim -pqm force.py -a
# manim -pqm rack_gear.py -a
# manim -pqm bevel_gear.py -a

PNG="_ManimCE_v0.17.2.png"

PRJROOT=".."
MEDIA=$PRJROOT/manim/media
SPUR=$PRJROOT/posts/involute-spur-gear
SNN=$PRJROOT/posts/introduction-spiking-nn
FL=$PRJROOT/posts/introduction-federated-learning
RL=$PRJROOT/posts/introduction-reinforcement-learning
BEVEL_1=$PRJROOT/posts/involute-bevel-gear-part-1
BEVEL_2=$PRJROOT/posts/involute-bevel-gear-part-2
BEVEL_3=$PRJROOT/posts/involute-bevel-gear-part-3

# force
force=$MEDIA/images/force/Gear2Gear$PNG
cp $force $SPUR/force.png

# involute
static_involute=$MEDIA/images/involute/Static$PNG
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
gearwithinterference=$MEDIA/images/rack_gear/GearWithInterference$PNG
gearwithoutinterference=$MEDIA/images/rack_gear/GearWithoutInterference$PNG
rackandgear=$MEDIA/images/rack_gear/RackAndGear$PNG
rack=$MEDIA/images/rack_gear/Rack$PNG
cp $gearwithinterference $SPUR/gear_with_interference.png
cp $gearwithoutinterference $SPUR/gear_without_interference.png
cp $rackandgear $SPUR/rack_and_gear.png
cp $rack $SPUR/rack.png

# signal
impulse=$MEDIA/videos/signal/720p30/Impulse.mp4
randomsignal=$MEDIA/images/signal/RandomSignal$PNG
cp $impulse $SNN
cp $randomsignal $SNN/random_signal.png

# federated_learning
aggregation=$MEDIA/videos/federated_learning/720p30/Aggregation.mp4
completetraining=$MEDIA/videos/federated_learning/720p30/CompleteTraining.mp4
downloadtoserver=$MEDIA/videos/federated_learning/720p30/DownloadToServer.mp4
restart=$MEDIA/videos/federated_learning/720p30/Restart.mp4
training=$MEDIA/videos/federated_learning/720p30/Training.mp4
uploadtocellphones=$MEDIA/videos/federated_learning/720p30/UploadToCellphones.mp4
selection=$MEDIA/images/federated_learning/Selection$PNG
serverandcellphones=$MEDIA/images/federated_learning/ServerAndCellphones$PNG
serverasnn=$MEDIA/images/federated_learning/ServerAsNN$PNG
cp $aggregation $FL
cp $completetraining $FL
cp $downloadtoserver $FL
cp $restart $FL
cp $training $FL
cp $uploadtocellphones $FL
cp $selection $FL/selection.png
cp $serverandcellphones $FL/server_and_cellphones.png
cp $serverasnn $FL/server_as_nn.png

# reinforcement_learning
context=$MEDIA/images/reinforcement_learning/Context$PNG
dice=$MEDIA/images/reinforcement_learning/Dice$PNG
expectedrewardinit=$MEDIA/images/reinforcement_learning/ExpectedRewardInit$PNG
simplified=$MEDIA/images/reinforcement_learning/Simplified$PNG
rlimage=$MEDIA/images/reinforcement_learning/Image$PNG
eat=$MEDIA/videos/reinforcement_learning/720p30/Eat.mp4
expectedrewardanimation=$MEDIA/videos/reinforcement_learning/720p30/ExpectedRewardAnimation.mp4
randompath=$MEDIA/videos/reinforcement_learning/720p30/RandomPath.mp4
rewardsdot=$MEDIA/videos/reinforcement_learning/720p30/RewardsDot.mp4
cp $context $RL/context.png
cp $dice $RL/dice.png
cp $expectedrewardinit $RL/expected_reward_init.png
cp $simplified $RL/simplified.png
cp $rlimage $RL/image.png
cp $eat $RL
cp $expectedrewardanimation $RL
cp $randompath $RL
cp $rewardsdot $RL

# New Manim version

PNG="_ManimCE_v0.19.0.png"

# bevel_gear
bevel_animation=$MEDIA/videos/bevel_gear/720p30/BevelAnimation.mp4
flat=$MEDIA/images/bevel_gear/BevelFlat$PNG
section=$MEDIA/images/bevel_gear/BevelGearSection$PNG
bevel_involute=$MEDIA/images/bevel_gear/BevelStaticInvolute$PNG
bevel_systems=$MEDIA/images/bevel_gear/SphericalRepr$PNG
bevel_params=$MEDIA/images/bevel_gear/SphericalRepr2$PNG
bevel_transmission=$MEDIA/images/bevel_gear/BevelTransmission$PNG
bevel_interference_params=$MEDIA/images/bevel_gear/BevelInterParams$PNG
# Part 1
cp $bevel_animation $BEVEL_1/involute-animation.mp4
cp $flat $BEVEL_1/flat.png
cp $section $BEVEL_1/section.png
cp $bevel_involute $BEVEL_1/involute.png
cp $bevel_transmission $BEVEL_1/transmission.png

# Part 2
cp $bevel_inter_params $BEVEL_2/interference-params.png
cp $bevel_transmission $BEVEL_2/transmission.png
cp $bevel_systems $BEVEL_2/systems.png
cp $bevel_params $BEVEL_2/params.png

# Part 3
cp $bevel_transmission $BEVEL_2/transmission.png
