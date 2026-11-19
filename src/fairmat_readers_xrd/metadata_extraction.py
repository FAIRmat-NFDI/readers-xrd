#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""Explicit, optional measurement metadata; independent of any archive schema."""

import hashlib
import re

import numpy as np
import pint

ureg = pint.get_application_registry()


def _text(value):
    return value.strip() if isinstance(value, str) and value.strip() else None


def _key(prefix, value):
    """Keep distinct source identities even when their display labels coincide."""
    label = re.sub(r'[^a-zA-Z0-9_]+', '_', str(value)).strip('_')[:40]
    digest = hashlib.sha256(str(value).encode()).hexdigest()[:10]
    return f'{prefix}_{label}_{digest}'


def _record(key, values, sources):
    values = {
        name: value
        for name, value in values.items()
        if value is not None and not (isinstance(value, str) and not value.strip())
    }
    return {
        'key': key,
        'values': values,
        'sources': {
            name: sources[name] if isinstance(sources[name], list) else [sources[name]]
            for name in values
        },
    }


def _same(left, right):
    if left is None or right is None:
        return left is right
    if hasattr(left, 'units'):
        if not hasattr(right, 'units') or not left.is_compatible_with(right):
            return False
        return np.array_equal(left.magnitude, right.to(left.units).magnitude)
    return np.array_equal(left, right)


def _collapse(records, logger=None):
    """Merge identical scan headers; omit conflicting fields instead of picking one."""
    groups = {}
    for record in records:
        groups.setdefault(record['key'], []).append(record)
    output = []
    for key, entries in groups.items():
        values, sources = {}, {}
        names = set().union(*(entry['values'] for entry in entries))
        for name in sorted(names):
            candidates = [entry['values'].get(name) for entry in entries]
            if not all(_same(candidates[0], value) for value in candidates[1:]):
                if logger is not None:
                    logger.warning(
                        f'Mixed scan metadata {key}/{name}; omitting common value.'
                    )
                continue
            values[name] = candidates[0]
            sources[name] = list(
                dict.fromkeys(
                    path for entry in entries for path in entry['sources'].get(name, [])
                )
            )
        if values:
            output.append({'key': key, 'values': values, 'sources': sources})
    return output


def _quantity(text, unit, path, logger, category=None):
    if text is None or not unit:
        return None
    try:
        value = float(text) * ureg(unit)
        if not np.isfinite(value.magnitude):
            raise ValueError('Non-finite value')
        if category and not value.is_compatible_with(category):
            raise ValueError('Incompatible unit')
        return value
    except (ValueError, TypeError, pint.PintError):
        if logger is not None:
            logger.warning(f'Ignoring invalid optional metadata at {path}.')
        return None


def xrdml_metadata(root, measurement, scans, ns, logger=None):
    details = {
        'users': [],
        'notes': [],
        'instrument_notes': [],
        'program': [],
        'detector': [],
    }
    sample_name = _text(root.findtext('xrdml:sample/xrdml:name', namespaces=ns))
    if sample_name:
        details['sample'] = [
            _record(
                'sample',
                {'sample_name': sample_name},
                {'sample_name': 'sample/name'},
            )
        ]
    for index, scan in enumerate(scans):
        base = f'xrdMeasurement/scan[{index + 1}]/header'
        author = _text(
            scan.findtext('xrdml:header/xrdml:author/xrdml:name', namespaces=ns)
        )
        if author:
            details['users'].append(
                _record(
                    _key('user', author),
                    {'name': author},
                    {'name': f'{base}/author/name'},
                )
            )
        program = scan.find('xrdml:header/xrdml:source/xrdml:applicationSoftware', ns)
        details['program'].append(
            _record(
                'acquisition_program',
                {
                    'name': _text(program.text) if program is not None else None,
                    'version': _text(program.get('version'))
                    if program is not None
                    else None,
                },
                {
                    'name': f'{base}/source/applicationSoftware',
                    'version': f'{base}/source/applicationSoftware/@version',
                },
            )
        )
        control = scan.find(
            'xrdml:header/xrdml:source/xrdml:instrumentControlSoftware', ns
        )
        if control is not None and _text(control.text):
            text = f'Instrument control software: {control.text.strip()}'
            if _text(control.get('version')):
                text += f' (version {control.get("version").strip()})'
            details['instrument_notes'].append(
                _record(
                    _key('control_software', text),
                    {'text': text},
                    {'text': f'{base}/source/instrumentControlSoftware'},
                )
            )
    # Empty headers must still participate in common-value checks.
    programs = details['program']
    if any(not record['values'].get('name') for record in programs):
        if any(record['values'] for record in programs) and logger is not None:
            logger.warning(
                'Incomplete scan acquisition program; omitting common metadata.'
            )
        details['program'] = []
    for parent, base in ((root, ''), (measurement, 'xrdMeasurement/')):
        for index, node in enumerate(parent.findall('xrdml:comment/xrdml:entry', ns)):
            text = _text(node.text)
            if text:
                details['notes'].append(
                    _record(
                        _key('comment', f'{base}{index}'),
                        {'text': text},
                        {'text': f'{base}comment/entry[{index + 1}]'},
                    )
                )
    aperture_tags = {
        'mask': 'width',
        'antiScatterSlit': 'height',
        'divergenceSlit': 'height',
        'receivingSlit': 'height',
    }
    for side in ('incidentBeamPath', 'diffractedBeamPath'):
        beam = measurement.find(f'xrdml:{side}', ns)
        if beam is None:
            continue
        occurrences = {}
        for node in beam:
            tag = node.tag.rsplit('}', 1)[-1]
            occurrences[tag] = occurrences.get(tag, 0) + 1
            index = occurrences[tag]
            base = f'xrdMeasurement/{side}/{tag}[{index}]'
            identity = node.get('id') or str(index)
            key = _key(f'{side}_{tag}', identity)
            name = _text(node.get('name'))
            if tag == 'detector':
                details['detector'].append(
                    _record(
                        'primary_detector',
                        {'description': name},
                        {'description': f'{base}/@name'},
                    )
                )
            elif tag in (
                'monochromator',
                'sollerSlit',
                *aperture_tags,
                'xRayMirror',
                'parallelPlateCollimator',
                'beamAttenuator',
            ):
                values, sources = {'name': name}, {'name': f'{base}/@name'}
                if tag == 'monochromator':
                    values['chemical_formula'] = _text(
                        node.findtext('xrdml:crystal', namespaces=ns)
                    )
                    sources['chemical_formula'] = f'{base}/crystal'
                    hkl = [
                        node.findtext(f'xrdml:hkl/xrdml:{axis}', namespaces=ns)
                        for axis in ('h', 'k', 'l')
                    ]
                    if all(value is not None for value in hkl):
                        try:
                            values['reflection'] = np.array(
                                [int(value) for value in hkl], dtype=np.int64
                            )
                            sources['reflection'] = [
                                f'{base}/hkl/{axis}' for axis in ('h', 'k', 'l')
                            ]
                        except ValueError:
                            if logger is not None:
                                logger.warning(
                                    f'Ignoring invalid reflection at {base}/hkl.'
                                )
                    family = 'monochromators'
                elif tag == 'sollerSlit' or tag in aperture_tags:
                    if tag == 'sollerSlit':
                        child_name = 'opening'
                        field = 'soller_angle'
                        unit = 'radian'
                        family = 'collimators'
                    else:
                        child_name = aperture_tags[tag]
                        field = 'size'
                        unit = 'meter'
                        family = 'apertures'
                    child = node.find(f'xrdml:{child_name}', ns)
                    values[field] = _quantity(
                        child.text if child is not None else None,
                        child.get('unit') if child is not None else None,
                        f'{base}/{child_name}',
                        logger,
                        unit,
                    )
                    sources[field] = f'{base}/{child_name}'
                    values['description'] = (
                        f'Reported {side}/{tag}; '
                        f'{child_name} is the supplied opening dimension.'
                    )
                    sources['description'] = base
                elif tag == 'beamAttenuator':
                    values['material'] = _text(
                        node.findtext('xrdml:material', namespaces=ns)
                    )
                    sources['material'] = f'{base}/material'
                    family = 'attenuators'
                else:
                    family = 'mirrors' if tag == 'xRayMirror' else 'collimators'
                details.setdefault(family, []).append(_record(key, values, sources))
    return {family: _collapse(records, logger) for family, records in details.items()}


def rasx_metadata(reader, logger=None):
    details = {
        'users': [],
        'notes': [],
        'instrument_notes': [],
        'detector': [],
        'sample': [],
        'positioners': [],
    }
    profiles = reader.meta[: len(reader.data)]
    for index, meta in enumerate(profiles):
        base = reader.metadata_paths[index]
        general = meta.get('GeneralInformation', {})
        operator = _text(general.get('Operator'))
        if operator:
            details['users'].append(
                _record(
                    _key('user', operator),
                    {'name': operator},
                    {'name': f'{base}/GeneralInformation/Operator'},
                )
            )
        for name in ('Comment', 'Memo', 'UserGroup', 'SystemName', 'Type', 'Version'):
            value = general.get(name)
            if value is not None and str(value).strip():
                text = str(value) if name in ('Comment', 'Memo') else f'{name}: {value}'
                details['notes'].append(
                    _record(
                        _key('rasx_note', f'{name}:{text}'),
                        {'text': text},
                        {'text': f'{base}/GeneralInformation/{name}'},
                    )
                )
        details['sample'].append(
            _record(
                'sample',
                {'sample_name': _text(general.get('SampleName'))},
                {'sample_name': f'{base}/GeneralInformation/SampleName'},
            )
        )
        hardware = meta.get('HardwareConfig', {})
        details['detector'].append(
            _record(
                'primary_detector',
                {'description': _text(hardware.get('Detector'))},
                {
                    'description': (
                        f'{base}/HWConfigurations/Categories/Detector/@SelectedUnit'
                    )
                },
            )
        )
        for name, value in hardware.get('optics_settings', {}).items():
            if _text(value):
                text = f'Optics {name}: {value.strip()}'
                details['instrument_notes'].append(
                    _record(
                        _key('rasx_optics', f'{name}:{value}'),
                        {'text': text},
                        {'text': f'{base}/HWConfigurations/Optics/{name}'},
                    )
                )
    # Only geometric, fixed, explicitly unit-bearing axes become positioner values.
    axes = {
        'ThetaS',
        'ThetaD',
        'TwoTheta',
        'Omega',
        'Chi',
        'Phi',
        'Z',
        'Zs',
        'Zr',
        'Ts',
        'Alpha',
        'Beta',
        'Rx',
        'Ry',
    }
    for name in sorted(axes):
        entries = [meta.get('Axes', {}).get(name) for meta in profiles]
        if not entries or any(
            axis is None or axis.State != 'Fixed' for axis in entries
        ):
            continue
        values = [
            _quantity(
                axis.Position,
                axis.Unit,
                f'{reader.metadata_paths[index]}/Axes/{name}/@Position',
                logger,
            )
            for index, axis in enumerate(entries)
        ]
        if any(
            value is None
            or not (
                value.is_compatible_with('radian') or value.is_compatible_with('meter')
            )
            for value in values
        ):
            continue
        unit = values[0].units
        if not all(value.is_compatible_with(unit) for value in values):
            continue
        details['positioners'].append(
            _record(
                _key('rasx_axis', name),
                {
                    'name': name,
                    'description': (
                        'Reported fixed axis position, in original profile order; '
                        'no target or time axis inferred.'
                    ),
                    'value': np.array([value.to(unit).magnitude for value in values])
                    * unit,
                },
                {
                    'name': (
                        f'{reader.metadata_paths[0]}/Axes/Axis[@Name="{name}"]/@Name'
                    ),
                    'description': (
                        f'{reader.metadata_paths[0]}/Axes/Axis[@Name="{name}"]/@State'
                    ),
                    'value': [
                        f'{reader.metadata_paths[i]}/Axes/'
                        f'Axis[@Name="{name}"]/@Position'
                        for i in range(len(entries))
                    ],
                },
            )
        )
    return {family: _collapse(records, logger) for family, records in details.items()}


def _as_list(value):
    return value if isinstance(value, list) else [value] if value else []


def brml_metadata(reader, logger=None):
    details = {'positioners': [], 'detector': []}
    for name, entries in reader.motors.items():
        if len(entries) != len(reader.data['ScanName']) or any(
            entry is None for entry in entries
        ):
            continue
        values = [
            _quantity(
                entry.get('Value'), entry.get('Unit'), f'Drives/{name}/Position', logger
            )
            for entry in entries
        ]
        if any(
            value is None
            or not (
                value.is_compatible_with('radian') or value.is_compatible_with('meter')
            )
            for value in values
        ):
            continue
        unit = values[0].units
        if not all(value.is_compatible_with(unit) for value in values):
            continue
        details['positioners'].append(
            _record(
                _key('brml_motor', name),
                {
                    'name': name,
                    'description': (
                        'Reported fixed drive position, in original RawData order.'
                    ),
                    'value': np.array([value.to(unit).magnitude for value in values])
                    * unit,
                },
                {
                    'name': (
                        f'FixedInformation/Drives/'
                        f'InfoData[@LogicName="{name}"]/@LogicName'
                    ),
                    'description': 'FixedInformation/Drives',
                    'value': [
                        f'{item["path"]}/FixedInformation/Drives/'
                        f'InfoData[@LogicName="{name}"]/Position/@Value'
                        for item in reader.fixed_information
                    ],
                },
            )
        )
    for item in reader.fixed_information:
        detector_count = len(details['detector'])
        for side in ('PrimaryTracks', 'SecondaryTracks'):
            tracks = item['value'].get('Instrument', {}).get(side, {})
            for track in _as_list(tracks.get('TrackInfoData')):
                for index, node in enumerate(
                    _as_list(track.get('MountedOptics', {}).get('InfoData'))
                ):
                    base = (
                        f'{item["path"]}/FixedInformation/Instrument/{side}/'
                        f'TrackInfoData[@LogicName="{track.get("@LogicName", "")}"]/'
                        f'MountedOptics/InfoData[{index + 1}]'
                    )
                    key = _key(
                        f'brml_{side}',
                        node.get('BeringInfo', {}).get('@ObjectPath')
                        or (
                            f'{track.get("@LogicName", "")}/'
                            f'{node.get("@LogicName", index)}'
                        ),
                    )
                    name = _text(node.get('@VisibleName'))
                    kind = node.get('@xsi:type', '')
                    class_path = node.get('BeringInfo', {}).get('@ClassPath', '')
                    values, sources = {'name': name}, {'name': f'{base}/@VisibleName'}
                    tube = node.get('Tube', {})
                    serial = _text(tube.get('SerialNumber', {}).get('@Value'))
                    if serial:
                        details.setdefault('source_fabrication', []).append(
                            _record(
                                'tube_fabrication',
                                {'serial_number': serial},
                                {'serial_number': f'{base}/Tube/SerialNumber/@Value'},
                            )
                        )
                    if '/Optic/Detector/' in class_path:
                        details['detector'].append(
                            _record(
                                'primary_detector',
                                {'description': name},
                                {'description': f'{base}/@VisibleName'},
                            )
                        )
                        continue
                    if kind in (
                        'FixedSlitInfoDataExt',
                        'VarSlitInfoData',
                        'DetectorSlitInfoData',
                    ):
                        field = 'SlitWidth' if 'SlitWidth' in node else 'Height'
                        size = node.get(field, {})
                        values['size'] = _quantity(
                            size.get('@Value'),
                            size.get('@Unit'),
                            f'{base}/{field}',
                            logger,
                            'meter',
                        )
                        sources['size'] = f'{base}/{field}/@Value'
                        values['description'] = (
                            f'Reported {side} aperture; opening dimension {field}.'
                        )
                        sources['description'] = base
                        family = 'apertures'
                    elif kind == 'GoebelMirrorInfoData':
                        family = 'mirrors'
                    else:
                        # Product names do not establish crystal material or geometry.
                        text = (
                            f'Mounted optic {side}/{node.get("@LogicName", index)}: '
                            f'{name or kind}; source type {kind}.'
                        )
                        details.setdefault('instrument_notes', []).append(
                            _record(key, {'text': text}, {'text': base})
                        )
                        continue
                    details.setdefault(family, []).append(_record(key, values, sources))
        if len(details['detector']) == detector_count:
            details['detector'].append(
                _record(
                    'primary_detector',
                    {'description': None},
                    {},
                )
            )
    return {family: _collapse(records, logger) for family, records in details.items()}


def raw_metadata(metadata):
    details = {'users': [], 'notes': [], 'program': []}
    user = _text(metadata.get('user'))
    if user:
        details['users'].append(
            _record(_key('user', user), {'name': user}, {'name': 'RAW4/USER'})
        )
    comment = _text(metadata.get('comment'))
    if comment:
        details['notes'].append(
            _record('raw4_comment', {'text': comment}, {'text': 'RAW4/COMMENT'})
        )
    creator = _text(metadata.get('creator'))
    if creator:
        details['program'].append(
            _record(
                'acquisition_program',
                {'name': creator, 'version': _text(metadata.get('creator_version'))},
                {'name': 'RAW4/CREATOR', 'version': 'RAW4/CREATOR_VERSION'},
            )
        )
    return details
