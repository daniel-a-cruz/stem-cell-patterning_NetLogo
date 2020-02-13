extensions [vid]
breed [ cells cell ]

globals [ altcount bhGain bhLoss blGain blLoss nhGain nhLoss ghGain ghLoss ]

turtles-own [ motion pluri divTrack diffTrack FGFR ERK GATA6 NANOG ]
patches-own [ FGF4 ]

to setup
  clear-all
  vid:reset-recorder

  set altcount 0
  set bhGain 0
  set bhLoss 0
  set blGain 0
  set blLoss 0
  set nhGain 0
  set nhLoss 0
  set ghGain 0
  set ghLoss 0

  set-default-shape turtles "circle"

  ask patches [
    ifelse (stochastic = "FGF4") or (stochastic = "All") [
      set FGF4 random (maxFGF4)
    ]
    [
      set FGF4 0
    ]
  ]

  create-cells (numNanogHigh) [
    setxy random-xcor random-ycor
    set motion true
    set pluri true
    set divTrack random (pluriMitosisThreshold)
    set diffTrack random (diffThresh)
    ifelse (stochastic = "FGFR+ERK") or (stochastic = "All") [
      set FGFR random 2
      set ERK random 2
    ]
    [
      set FGFR 0
      set ERK 0
    ]
    ; The setup below sets this cell as Gata6 Low
    set GATA6 0
    set NANOG 1
    set color lime
  ]

  create-cells (numGataHigh) [
    setxy random-xcor random-ycor
    set motion true
    set pluri true
    set divTrack random (pluriMitosisThreshold)
    set diffTrack random (diffThresh)
    ifelse (stochastic = "FGFR+ERK") or (stochastic = "All") [
      set FGFR random 2
      set ERK random 2
    ]
    [
      set FGFR 0
      set ERK 0
    ]
    ; The setup below sets this cell as Gata6 High
    set GATA6 1
    set NANOG 0
    set color white
  ]

  reset-ticks
end

to go
  set bhGain 0
  set bhLoss 0
  set blGain 0
  set blLoss 0
  set nhGain 0
  set nhLoss 0
  set ghGain 0
  set ghLoss 0

  ; Synchronized updates, currently not in use
  ;ifelse synchronization [
  ;  ask cells [ set divTrack divTrack + 5 ]
  ;  ask cells [ single-cell-move ]
  ;  ask cells [
  ;    ifelse diffInteract [
  ;      FE-pathway
  ;    ]
  ;    [
  ;      if pluri [ FE-pathway ]
  ;    ]
  ;  ]
  ;  ask cells [ diff-low-surrounded ]
  ;  ask cells [ cell-division ]
  ;  ask cells [ color-update ]
  ;]

  ask cells [
    set divTrack divTrack + 5
    single-cell-move
    ifelse diffInteract [
      FE-pathway
    ]
    [
      if pluri [ FE-pathway ]
    ]
    diff-low-surrounded
    ;spontaneous-diff
    cell-division
    color-update
    repel
  ]

  ask patches [
    if FGF4 > 0 [
      set FGF4 FGF4 - 1
    ]
  ]

  tick

  let ylw count cells with [ color = yellow ]
  let blu count cells with [ color = sky ]
  set altcount (ylw + blu)
end

; Cone movement, currently not in use
to move [ dist ]
  if dist >= .25 [
    let obs count cells in-cone 1 60
    ifelse obs <= crowdThreshold [
      fd dist
    ]
    [
      let angle random-float 30
      let temp random 2
      ifelse temp = 1 [
        rt 30 + angle
      ]
      [
        lt 30 + angle
      ]
      move dist / 2
    ]
  ]
end

to single-cell-move
  if motion [
    ;ifelse (NANOG = 0) and (GATA6 = 1) and pluri [
    ifelse (GATA6 = 1) and pluri [
      ifelse guyeMove [
        let nearestDiff min-one-of cells with [pluri = false ] [ distance myself ]
        ifelse nearestDiff != nobody [
          face nearestDiff
          fd cellSpeed
        ]
        [
          rt random-float 360
          fd cellSpeed
        ]
      ]
      [
        rt random-float 360
        fd cellSpeed
      ]
    ]
    [
      ifelse pluri and nClustering [
        let nearestLow min-one-of cells with [ NANOG = 0 ] [ distance myself ]
        ifelse nearestLow != nobody [
          face nearestLow
          fd cellSpeed
        ]
        [
          rt random-float 360
          fd cellSpeed
        ]
      ]
      [
        rt random-float 360
        fd cellSpeed
      ]
    ]

    if not pluri [
      if any? other cells with [ pluri = false ] in-radius 1.0 [
        set motion false
      ]
    ]

    if (NANOG = 1) and (GATA6 = 0) and pluri [
      if any? other cells with [ (NANOG = 1) and (GATA6 = 0) and (pluri = true) ] in-radius 1.0 [
        set motion false
      ]
    ]
  ]
end

to differentiate
  set motion true
  set pluri false
  set color red
  set GATA6 1
  set NANOG 0
end

to FE-pathway
  let boolFGF4 0
  let tempFGFR FGFR
  let tempERK ERK
  let tempNANOG NANOG
  let tempGATA6 GATA6

  ask patch-here [
    if FGF4 > 0 [
      set boolFGF4 1
    ]
    ; Boolean function for FGF4
    if (tempNANOG = 1) and (FGF4 < maxFGF4) [
      set FGF4 FGF4 + 1
    ]
  ]

  ifelse wDiagram = "FGFR->ERK" [
    ; Boolean function for ERK depends only on FGFR
    set ERK tempFGFR
    ; Boolean function for FGFR based on FGF4 and GATA6
    ifelse func = "AND" [
      ; AND
      set FGFR (boolFGF4 * tempGATA6)
    ]
    ; OR
    [
      set FGFR ((boolFGF4 + tempGATA6) + (boolFGF4 * tempGATA6)) mod 2
    ]
  ]
  [
    ; Boolean function for FGFR depends only on GATA6
    set FGFR tempGATA6
    ; Boolean function for ERK based on FGF4 and FGFR
    ifelse func = "AND" [
      ; AND
      set ERK (boolFGF4 * tempFGFR)
    ]
    ; OR
    [
      set ERK ((boolFGF4 + tempFGFR) + (boolFGF4 * tempFGFR)) mod 2
    ]
  ]


  ; Other Boolean functions
  set GATA6 (tempNANOG + 1) mod 2
  set NANOG ((tempERK + 1) mod 2) * ((tempGATA6 + 1) mod 2)

  if (tempFGFR = 0) and (FGFR = 1) [
    ask patch-here [ set FGF4 FGF4 - 1 ]
  ]

  ifelse (GATA6 = 1) and pluri [
    set diffTrack diffTrack + 1
    if diffTrack >= diffThresh [ differentiate ]
  ]
  [
    set diffTrack 0
  ]

  ; Added to prevent cells that are not Nanog high & Gata6 low
  ; from clustering if there was a change
  if (tempNANOG = 0) or (tempGATA6 = 1) and pluri [
    set motion true
  ]

  ; Nanog High tracking

  if (NANOG = 1) and (GATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set nhGain nhGain + 1
  ]

  if (tempNANOG = 1) and (tempGATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set nhLoss nhLoss + 1
  ]

  ; Gata6 High tracking

  if (NANOG = 0) and (GATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set ghGain ghGain + 1
  ]

  if (tempNANOG = 0) and (tempGATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set ghLoss ghLoss + 1
  ]

  ; Both Bigh tracking

  if (NANOG = 1) and (GATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set bhGain bhGain + 1
  ]

  if (tempNANOG = 1) and (tempGATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set bhLoss bhLoss + 1
  ]


  ; Both Low tracking


  if (NANOG = 0) and (GATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set blGain blGain + 1
  ]

  if (tempNANOG = 0) and (tempGATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set blLoss blLoss + 1
  ]

end

to diff-low-surrounded
  if diffLowSurr [
    ;if (NANOG = 1) and (GATA6 = 0) and pluri [
    if (GATA6 = 0) and pluri [
      let crowd count cells with [ pluri = false ] in-radius interactionDistance
      if crowd > crowdThreshold [
        set diffTrack diffTrack + 1

        ;if diffTrack >= diffThresh [ differentiate ]

        ; Spontaneous differentiation
        differentiate
      ]
    ]
  ]
end

; Spontaneous differentiation, currently not in use
to spontaneous-diff
;  if sponDiff != "Off" and pluri [
;    let prob random-float 100
;    if sponDiff = "Gata6 Low + Nanog High" and (NANOG = 1) and (GATA6 = 0) and (prob < sponProb) [
;      differentiate
;    ]

;    if sponDiff = "Gata6 Low"  and (GATA6 = 0) and (prob < sponProb) [
;      differentiate
;    ]

;    if sponDiff = "All" and (prob < sponProb) [
;      differentiate
;    ]
;  ]
end

to cell-division
  if not motion [
    ifelse not pluri [
      if (divTrack > diffMitosisThreshold) [
        set divTrack (divTrack / 2)
        hatch 1 [
          rt random-float 360
          fd 1
          repel
        ]
      ]
    ]
    [
      if (divTrack > pluriMitosisThreshold) [
        set divTrack (divTrack / 2)
        hatch 1 [
          rt random-float 360
          fd 1
          repel
        ]
      ]
    ]
  ]
end

to color-update
  ifelse not pluri [
    set color red
  ]
  [
    if (NANOG = 1) and (GATA6 = 0) [
      set color lime
    ]
    if (NANOG = 0) and (GATA6 = 1) [
      set color white
    ]
    if (NANOG = 1) and (GATA6 = 1) [
      set color yellow
    ]
    if (NANOG = 0) and (GATA6 = 0) [
      set color sky
    ]
  ]
end

to repel
  let too-near one-of other cells in-radius 0.9
  if too-near != nobody [
    face too-near
    fd -1 * (cellSpeed / 5)
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
275
10
889
625
-1
-1
6.0
1
15
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
steps
2.0

SLIDER
5
135
265
168
numGataHigh
numGataHigh
0
1000
500.0
10
1
NIL
HORIZONTAL

SLIDER
5
205
265
238
numNanogHigh
numNanogHigh
0
1000
500.0
10
1
NIL
HORIZONTAL

BUTTON
900
680
973
713
Initialize
setup\n\nif record-type = \"View\" [\nvid:start-recorder\nvid:record-view\n]\n\nif record-type = \"Interface\" [\nvid:start-recorder\nvid:record-interface\n]    
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1065
680
1150
713
Run
go\n\nif record-type = \"View\" [\nvid:record-view\n]\n\nif record-type = \"Interface\" [\nvid:record-interface\n] 
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
5
420
265
453
interactionDistance
interactionDistance
1
5
2.0
1
1
NIL
HORIZONTAL

MONITOR
700
635
780
680
Differentiated
count cells with [ color = red ]
17
1
11

SLIDER
5
275
265
308
pluriMitosisThreshold
pluriMitosisThreshold
0
500
75.0
5
1
NIL
HORIZONTAL

SLIDER
5
345
265
378
diffMitosisThreshold
diffMitosisThreshold
0
500
100.0
5
1
NIL
HORIZONTAL

MONITOR
280
635
360
680
Nanog High
count cells with [ color = lime ]
17
1
11

BUTTON
982
680
1057
713
Step
go\n\nif record-type = \"View\" [\nvid:record-view\n]\n\nif record-type = \"Interface\" [\nvid:record-interface\n] 
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1045
725
1150
765
Save Video
if record-type = \"View\" or record-type = \"Interface\" [\nvid:save-recording user-new-file\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
5
510
265
543
crowdThreshold
crowdThreshold
5
30
15.0
5
1
NIL
HORIZONTAL

CHOOSER
900
725
1035
770
record-type
record-type
"None" "View" "Interface"
1

MONITOR
385
635
465
680
Gata6 High
count cells with [ color = white ]
17
1
11

TEXTBOX
75
10
180
40
Parameters
20
0.0
1

SLIDER
5
60
265
93
cellSpeed
cellSpeed
0
5
2.0
.5
1
NIL
HORIZONTAL

TEXTBOX
965
10
1080
40
Assumptions
20
0.0
1

SWITCH
900
245
1150
278
guyeMove
guyeMove
1
1
-1000

TEXTBOX
905
50
1150
68
Stochastic FGFR+ERK and/or FGF4?
13
0.0
1

TEXTBOX
905
220
1150
238
Movement: Random or Guye model?
13
0.0
1

TEXTBOX
905
460
1150
491
Pluripotent Gata6 Low cells differentiate within differentiated subpopulation?
13
0.0
1

SWITCH
900
500
1150
533
diffLowSurr
diffLowSurr
1
1
-1000

TEXTBOX
905
550
1150
581
Do differentiated cells interact with the FGF/ERK pathway (via FGF4)?
13
0.0
1

SWITCH
900
590
1150
623
diffInteract
diffInteract
1
1
-1000

TEXTBOX
985
640
1060
660
Controls
20
0.0
1

TEXTBOX
10
40
260
58
Cell speed
13
0.0
1

TEXTBOX
10
115
260
133
Number of initial Gata6 High Cells
13
0.0
1

TEXTBOX
10
185
260
205
Number of intiial Nanog High Cells
13
0.0
1

TEXTBOX
10
325
260
343
Mitosis threshold for pluripotent cells
13
0.0
1

TEXTBOX
10
255
260
273
Mitosis threshold for differentiated cells\t
13
0.0
1

TEXTBOX
10
400
265
418
Interaction distance for determining a crowd
13
0.0
1

TEXTBOX
10
470
260
500
Crowd threshold when differentiating Gata6 Low cells near differentiated cells
13
0.0
1

TEXTBOX
10
565
250
596
Threshold (in steps) for differentiating cells
13
0.0
1

SLIDER
5
585
265
618
diffThresh
diffThresh
0
30
30.0
3
1
NIL
HORIZONTAL

MONITOR
490
635
570
680
Both High
count cells with [ color = yellow ]
17
1
11

TEXTBOX
1270
10
1355
31
Monitors
20
0.0
1

SLIDER
5
660
265
693
maxFGF4
maxFGF4
1
5
5.0
1
1
NIL
HORIZONTAL

TEXTBOX
10
640
260
658
Maximum FGF4 stored on a patch
13
0.0
1

MONITOR
595
635
677
680
Both Low
count cells with [ color = sky ]
17
1
11

CHOOSER
900
75
1145
120
stochastic
stochastic
"Off" "FGFR+ERK" "FGF4" "All"
0

CHOOSER
900
160
1145
205
wDiagram
wDiagram
"FGFR->ERK" "GATA6->FGFR"
0

TEXTBOX
905
135
1140
153
Spontaneous differentation?
13
0.0
1

TEXTBOX
285
685
355
705
Green
15
0.0
1

TEXTBOX
390
685
460
705
White
15
0.0
1

TEXTBOX
495
685
565
705
Yellow
15
0.0
1

TEXTBOX
600
685
670
703
Blue
15
0.0
1

TEXTBOX
705
685
775
703
Red
15
0.0
1

MONITOR
1305
130
1395
175
Gata6 High Loss
ghLoss
17
1
11

MONITOR
1305
75
1395
120
Nanog High Loss
nhLoss
17
1
11

MONITOR
1195
75
1285
120
Nanog High Gain
nhGain
17
1
11

MONITOR
1195
130
1285
175
Gata6 High Gain
ghGain
17
1
11

MONITOR
805
635
885
680
Alt. Cells
altcount
17
1
11

TEXTBOX
810
685
880
721
Blue + Yellow
15
0.0
1

MONITOR
1195
185
1285
230
Both High Gain
bhGain
17
1
11

MONITOR
1305
185
1395
230
Both High Loss
bhLoss
17
1
11

MONITOR
1195
240
1285
285
Both Low Gain
blGain
17
1
11

MONITOR
1305
240
1395
285
Both Low Loss
blLoss
17
1
11

TEXTBOX
905
370
1145
388
FGFR/ERK function: 'AND' or 'OR'?
13
0.0
1

CHOOSER
900
395
1150
440
func
func
"AND" "OR"
1

TEXTBOX
905
295
1140
313
Nanog High clustering?
13
0.0
1

SWITCH
900
320
1150
353
nClustering
nClustering
1
1
-1000

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
