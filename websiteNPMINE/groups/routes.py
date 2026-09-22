from flask import Blueprint, abort, flash, redirect, render_template, url_for
from flask_login import current_user, login_required
from sqlalchemy import func, or_

from websiteNPMINE import db
from websiteNPMINE.groups.forms import CompoundGroupForm, GroupForm, InviteMemberForm
from websiteNPMINE.models import AccountGroup, Accounts, Compounds, Group

groups = Blueprint('groups', __name__)


@groups.route('/groups', methods=['GET', 'POST'])
@login_required
def index():
    group_form = GroupForm()
    invite_form = InviteMemberForm()

    if group_form.validate_on_submit():
        group = Group(
            name=group_form.name.data.strip(),
            owner=current_user
        )
        db.session.add(group)
        db.session.flush()
        db.session.add(AccountGroup(group=group, account=current_user, role='editor'))
        db.session.commit()

        flash(f"Group '{group.name}' created.", 'success')
        return redirect(url_for('groups.index'))
    elif group_form.is_submitted():
        for errors in group_form.errors.values():
            for error in errors:
                flash(error, 'error')

    user_groups = (
        Group.query
        .filter(or_(Group.user_id == current_user.id, Group.members.any(Accounts.id == current_user.id)))
        .order_by(Group.name.asc())
        .all()
    )

    return render_template(
        'groups.html',
        logged_in=current_user.is_authenticated,
        group_form=group_form,
        invite_form=invite_form,
        groups=user_groups
    )


@groups.route('/groups/<int:group_id>/members', methods=['POST'])
@login_required
def add_member(group_id):
    group = Group.query.get_or_404(group_id)
    if group.user_id != current_user.id:
        abort(403)

    form = InviteMemberForm()
    if not form.validate_on_submit():
        for errors in form.errors.values():
            for error in errors:
                flash(error, 'error')
        return redirect(url_for('groups.index'))

    email = form.email.data.strip().lower()
    account = Accounts.query.filter(func.lower(Accounts.email) == email).first()

    if account is None:
        flash('No account exists with that email.', 'error')
        return redirect(url_for('groups.index'))

    existing_membership = AccountGroup.query.filter_by(group_id=group.id, account_id=account.id).first()
    if existing_membership:
        flash(f'{account.email} is already in this group.', 'info')
        return redirect(url_for('groups.index'))

    db.session.add(AccountGroup(group=group, account=account, role=form.role.data))
    db.session.commit()

    flash(f'{account.email} was added to {group.name} as {form.role.data}.', 'success')
    return redirect(url_for('groups.index'))


@groups.route('/compounds/<int:compound_id>/groups', methods=['POST'])
@login_required
def add_compound(compound_id):
    compound = Compounds.active().filter_by(id=compound_id).first_or_404()
    if current_user.role_id != 1 and compound.user_id != current_user.id:
        abort(403)

    linked_group_ids = {group.id for group in compound.groups.all()}
    available_groups_query = (
        Group.query
        .filter(or_(Group.user_id == current_user.id, Group.members.any(Accounts.id == current_user.id)))
    )
    if linked_group_ids:
        available_groups_query = available_groups_query.filter(~Group.id.in_(linked_group_ids))

    available_groups = available_groups_query.order_by(Group.name.asc()).all()

    form = CompoundGroupForm()
    form.group_id.choices = [(group.id, group.name) for group in available_groups]

    if not form.validate_on_submit():
        for errors in form.errors.values():
            for error in errors:
                flash(error, 'error')
        return redirect(url_for('main.compound', compound_id=compound.id))

    group = next((group for group in available_groups if group.id == form.group_id.data), None)
    if group is None:
        flash('Select a valid group.', 'error')
        return redirect(url_for('main.compound', compound_id=compound.id))

    compound.groups.append(group)
    db.session.commit()

    flash(f"'{compound.compound_name}' was added to {group.name}.", 'success')
    return redirect(url_for('main.compound', compound_id=compound.id))
