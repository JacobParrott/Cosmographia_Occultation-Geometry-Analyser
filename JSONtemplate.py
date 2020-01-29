#FORM A TEMPLATE JSON FILE
# Profile begining
JSONstart = '''
{
  "version": "1.0",
  "name": "Occultation Profile",
  "items": [
'''
# empty string required to grow from
JSONiterated = ''

# this one has the comma on it. lenght of profile adjust the iteration count of this
JSONiterable = '''

{
  "name": "Profile Iterate",
  "center": "Mars",
  "trajectory": {
    "type": "FixedPoint",
    "position": [
      0.0,
      0.0,
      0.0
    ]
  },
  "bodyFrame": {
    "type": "BodyFixed",
    "body": "Mars"
  },
  "geometry": {
    "type": "ParticleSystem",
    "emitters": [
      {
        "texture": "gaussian.jpg",
        "generator": {
          "type": "Strip",
          "states": []
        },
        "velocityVariation": 8,
        "trace": 2,
        "spawnRate": 20,
        "lifetime": 30,
        "startSize": 5,
        "endSize": 10,

        "colors": [
          "#00ff00",
          0.000,
          "#ff0000",
          0.00,
          "#0000ff",
          0.005,
          "#fffff",
          0.005
        ],
        "startTime":"2020-01-01 12:00:00",
        "endTime": "2020-01-01 14:30:00"
      }
    ]
  }
},

'''
# must be placed at the end of the json stack[ doesnt have a comma for the last item]
JSONend = '''

{
  "name": "Final Profile [high in the ionosphere]",
  "center": "Mars",
  "trajectory": {
    "type": "FixedPoint",
    "position": [
      0.0,
      0.0,
      0.0
    ]
  },
  "bodyFrame": {
    "type": "BodyFixed",
    "body": "Mars"
  },
  "geometry": {
    "type": "ParticleSystem",
    "emitters": [
      {
        "texture": "gaussian.jpg",
        "generator": {
          "type": "Strip",
          "states": []
        },
        "velocityVariation": 2,
        "trace": 2,
        "spawnRate": 100,
        "lifetime": 30,
        "startSize": 5,
        "endSize": 100,

        "colors": [
          "#111111",
          0.0,
          "#ffffff",
          0.00,
          "#111111",
          0.045,
          "#0000ff",
          0.0
        ],
        "startTime":"2020-01-01 12:00:00",
        "endTime": "2020-01-01 14:30:00"
      }
    ]
  }
}
]
}

'''