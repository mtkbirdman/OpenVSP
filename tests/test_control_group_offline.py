"""Offline regression tests: group enablement cannot be inferred from a set angle."""
import pytest
from openvsppy.core.util import set_control_surface


class FakeVSP:
    def __init__(self, *, assigned=True, active=False):
        self.names = ['AILERON_GROUP']
        self.surfaces = (['WingGeom_Surf0_AILERON', 'WingGeom_Surf1_AILERON']
                         if assigned else [])
        self.available = ([] if assigned else
                          ['WingGeom_Surf0_AILERON', 'WingGeom_Surf1_AILERON'])
        self.parms = {'ActiveFlag': float(active), 'DeflectionAngle': 0.,
                      'Surf_SUB0_0_Gain': 0., 'Surf_SUB0_1_Gain': 0.}
        self.calls_to_add = []

    def FindGeomsWithName(self, name):
        return ['GEOM0'] if name == 'WingGeom' else []

    def GetSubSurf(self, geom, sub_id):
        return 'SUB0' if geom == 'GEOM0' and sub_id == 0 else ''

    def GetSubSurfName(self, sub_id):
        return 'AILERON'

    def GetNumControlSurfaceGroups(self):
        return len(self.names)

    def GetVSPAEROControlGroupName(self, idx):
        return self.names[idx]

    def CreateVSPAEROControlSurfaceGroup(self):
        self.names.append('')
        return len(self.names) - 1

    def SetVSPAEROControlGroupName(self, name, idx):
        self.names[idx] = name

    def GetActiveCSNameVec(self, idx):
        return list(self.surfaces)

    def GetAvailableCSNameVec(self, idx):
        return list(self.available)

    def AddSelectedToCSGroup(self, selected, idx):
        self.calls_to_add.append(selected)
        assert all(1 <= i <= len(self.available) for i in selected)
        self.surfaces.extend(self.available[i - 1] for i in selected)
        for i in sorted(selected, reverse=True):
            del self.available[i - 1]

    def FindContainer(self, name, idx):
        return 'SETTINGS' if name == 'VSPAEROSettings' else ''

    def FindParm(self, settings, name, group):
        return name if settings == 'SETTINGS' and name in self.parms and group == 'ControlSurfaceGroup_0' else ''

    def SetParmVal(self, param, value):
        self.parms[param] = value

    def GetParmVal(self, param):
        return self.parms[param]

    def Update(self):
        pass


@pytest.mark.parametrize('assigned', [True, False])
def test_set_control_surface_enables_group_and_preserves_membership(assigned):
    vsp = FakeVSP(assigned=assigned)
    for angle in (10., -10., 0.):
        parm_id = set_control_surface(vsp, 'WingGeom', angle, 'AILERON_GROUP',
                                      gains=(1., 1.), verbose=0)
        assert parm_id == 'DeflectionAngle'
        assert vsp.parms['ActiveFlag'] == 1
        assert vsp.parms['DeflectionAngle'] == angle
        assert vsp.parms['Surf_SUB0_0_Gain'] == 1
        assert vsp.parms['Surf_SUB0_1_Gain'] == 1
    assert len(vsp.calls_to_add) == (0 if assigned else 1)
    assert len(vsp.surfaces) == 2


def test_set_control_surface_rejects_unavailable_control_surfaces():
    vsp = FakeVSP(assigned=False)
    vsp.available.clear()
    with pytest.raises(RuntimeError, match='expected 2 active'):
        set_control_surface(vsp, 'WingGeom', 10., 'AILERON_GROUP', gains=(1, 1), verbose=0)
